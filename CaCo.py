import os
import sys
import json
import lzma
import argparse
import tempfile
import shutil
import subprocess          
import multiprocessing as mp
from glob import glob
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm
import pyrodigal


# ---------- Helper: get file path relative to script ----------
def get_file_path(infile):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, infile)


# ---------- Download dbCAN HMM database ----------
def download_db():
    ofolder = get_file_path('data/')
    os.makedirs(ofolder, exist_ok=True)
    file_path = f'{ofolder}/dbcan.hmm'
    result = subprocess.run(
        ['wget', '-O', file_path,
         'https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-V11.txt'],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return f"Download failed: {result.stderr}"
    # Validate header/footer
    try:
        with open(file_path, 'r') as f:
            first = f.readline().strip()
            if not first.startswith("HMMER3/b [3.0 | March 2010]"):
                return "Database header incorrect"
            f.seek(0, os.SEEK_END)
            size = f.tell()
            if size < 2:
                return "File too short"
            f.seek(max(0, size - 10), os.SEEK_SET)
            tail = f.read().rstrip()
            if not tail.endswith("//"):
                return "Missing end marker"
    except Exception as e:
        return f"Error reading file: {e}"
    return "Download successful"


# ---------- Gene prediction with Pyrodigal (parallel) ----------
def predict_genes_pyrodigal(infile, outdir):
    """Predict genes from nucleotide genome, return path to .faa."""
    outfile = os.path.join(outdir, os.path.basename(infile).rsplit('.', 1)[0] + '.faa')
    if os.path.exists(outfile):
        return outfile
    with open(infile, 'rb') as f:
        seq = f.read()
    # Single genome mode (set meta=True for metagenomes)
    genes = pyrodigal.GeneFinder(meta=False)
    genes.train(seq)
    with open(outfile, 'w') as out:
        for idx, pred in enumerate(genes.find_genes(seq), 1):
            prot = pred.translate()
            out.write(f">{os.path.basename(infile)}_gene_{idx}\n{prot}\n")
    return outfile


# ---------- HMM search with hmmsearch (parallel) ----------
def run_hmmsearch(genome_faa, db, outdir, cpu_per_task=2):
    outfile = os.path.join(outdir, os.path.basename(genome_faa).replace('.faa', '.dbcan'))
    if os.path.exists(outfile):
        return outfile
    cmd = [
        'hmmsearch',
        '--cpu', str(cpu_per_task),
        '-o', '/dev/null',
        '--domtblout', outfile,
        db, genome_faa
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return outfile


def parallel_hmmsearch(faa_files, db, outdir, num_workers, cpu_per_task=2):
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = {ex.submit(run_hmmsearch, faa, db, outdir, cpu_per_task): faa
                   for faa in faa_files}
        results = []
        for fut in tqdm(as_completed(futures), total=len(futures), desc="HMM search"):
            results.append(fut.result())
    return results


# ---------- Parse hmmsearch output (parallel) ----------
def parse_hmmsearch_output(infile, outdir):
    parsed_file = f'{infile}.parsed'
    if os.path.exists(parsed_file):
        return parsed_file

    min_score = 25
    min_coverage = 0.35
    max_evalue = 1e-15
    best_hits = {}
    current_query = None
    current_best = None

    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 18:
                continue
            qname = fields[0]
            hname = fields[3].split('.')[0]
            evalue = float(fields[6])
            score = float(fields[7])
            hstart = int(fields[15])
            hend = int(fields[16])
            qlen = int(fields[2])
            cov = min(1.0, (hend - hstart + 1) / qlen)

            if qname != current_query:
                if current_query is not None and current_best is not None:
                    if (current_best[5] >= min_coverage and
                        current_best[1] <= max_evalue and
                        current_best[2] >= min_score):
                        best_hits[current_query] = current_best
                current_query = qname
                current_best = None

            if (current_best is None or score > current_best[2] or
                (score == current_best[2] and evalue < current_best[1]) or
                (score == current_best[2] and evalue == current_best[1] and cov > current_best[5])):
                current_best = (hname, evalue, score, hstart, hend, cov)

    if current_query is not None and current_best is not None:
        if (current_best[5] >= min_coverage and
            current_best[1] <= max_evalue and
            current_best[2] >= min_score):
            best_hits[current_query] = current_best

    df = pd.DataFrame.from_dict(best_hits, orient='index',
                                columns=['Hit_Name','Hit_Evalue','Hit_Score',
                                         'Hit_Start','Hit_End','Hit_Coverage'])
    df.to_csv(parsed_file, sep='\t', index_label='Query_Domain')
    return parsed_file


def parallel_parse(dbcan_files, outdir, num_workers):
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = {ex.submit(parse_hmmsearch_output, f, outdir): f for f in dbcan_files}
        results = []
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Parsing HMM"):
            results.append(fut.result())
    return results


# ---------- Feature extraction (parallel) ----------
def extract_features_one(parsed_file, subs_dict):
    name = os.path.basename(parsed_file).replace('.dbcan.parsed', '')
    df = pd.read_table(parsed_file)
    families = set(df.Hit_Name)
    families = {x.split('_')[0] for x in families}
    sorted_fam = sorted(families)
    subs_set = set()
    for fam in families:
        if fam in subs_dict:
            subs_set.update(subs_dict[fam].split(', '))
    sorted_subs = sorted(subs_set)
    return name, sorted_fam, sorted_subs


def extract_features_parallel(parsed_files, subs_dict, num_workers):
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = {ex.submit(extract_features_one, pf, subs_dict): pf for pf in parsed_files}
        fam_records = []
        sub_records = []
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Extracting features"):
            name, fams, subs = fut.result()
            fam_records.append({'genome': name, 'families': ', '.join(fams)})
            sub_records.append({'genome': name, 'substrates': ', '.join(subs)})
    return pd.DataFrame(fam_records), pd.DataFrame(sub_records)


# ---------- Probability overlap ----------
def probability_overlap(n1, n2, k, trials=1000):
    import random
    random.seed(42)
    M = 60
    cnt = 0
    for _ in range(trials):
        s1 = set(random.sample(range(M), min(n1, M)))
        s2 = set(random.sample(range(M), min(n2, M)))
        if len(s1 & s2) >= k:
            cnt += 1
    return cnt / trials


def probcomp(m, n, k):
    if (m+n) <= k or (m==0 and n==0):
        return 0, 0, 1.0
    if m==0 or n==0:
        return 0, 0, 1.0
    if not (m >= k and n >= k):
        return 'n.a.', 'n.a.', 'n.a.'
    maxcomp = min(m,n) / (m+n-min(m,n))
    comp = k/(m+n-k)
    relcomp = comp / maxcomp if maxcomp != 0 else 0
    prob = probability_overlap(m, n, k)
    return comp, relcomp, prob


# ---------- Pairwise RPS (parallel with chunking) ----------
def compute_rps_chunk(pair_chunk, df_sub, temp_dir):
    """Process a list of pairs and write results to a temp file."""
    temp_out = os.path.join(temp_dir, f"rps_{os.getpid()}.tsv")
    with open(temp_out, 'w') as f:
        f.write("genome1\tgenome2\tset1\tset2\tintersection\tcompetition\trelcomp\tprob\tRPS\trelRPS\n")
        for g1, g2 in pair_chunk:
            row1 = df_sub[df_sub.genome == g1]
            row2 = df_sub[df_sub.genome == g2]
            if row1.empty or row2.empty:
                continue
            set1 = set(row1.iloc[0]['substrates'].split(', '))
            set2 = set(row2.iloc[0]['substrates'].split(', '))
            intersect = len(set1 & set2)
            m, n = len(set1), len(set2)
            competition, relcomp, prob = probcomp(m, n, intersect)
            if isinstance(competition, str):
                continue
            RPS = 1 - 2 * competition
            relRPS = 1 - 2 * relcomp
            f.write(f"{g1}\t{g2}\t{m}\t{n}\t{intersect}\t{competition}\t{relcomp}\t{prob}\t{RPS}\t{relRPS}\n")
    return temp_out


def compute_rps_parallel(df_sub, num_workers, output_file):
    genomes = df_sub['genome'].tolist()
    pairs = list(combinations(genomes, 2))
    if not pairs:
        return

    chunk_size = max(1, len(pairs) // num_workers)
    chunks = [pairs[i:i+chunk_size] for i in range(0, len(pairs), chunk_size)]

    temp_dir = tempfile.mkdtemp(prefix="rps_chunks_")
    try:
        with ProcessPoolExecutor(max_workers=num_workers) as ex:
            futures = [ex.submit(compute_rps_chunk, ch, df_sub, temp_dir) for ch in chunks]
            temp_files = []
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Computing RPS"):
                temp_files.append(fut.result())

        # Merge all temp files into final compressed output
        with lzma.open(output_file, 'wt') as out_f:
            header_written = False
            for tf in temp_files:
                with open(tf, 'r') as inf:
                    if not header_written:
                        out_f.write(inf.readline())
                        header_written = True
                    else:
                        inf.readline()  # skip header
                    for line in inf:
                        out_f.write(line)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# ---------- Main orchestration ----------
def main(mode, genomes_list, temp_dir, output_file, db, subs_dict, num_workers):
    os.makedirs(temp_dir, exist_ok=True)

    # 1. Gene prediction (if nucleotides)
    if mode == 'from_nucleotides':
        print("Gene prediction with Pyrodigal (parallel)...")
        faa_files = []
        with ProcessPoolExecutor(max_workers=num_workers) as ex:
            futures = {ex.submit(predict_genes_pyrodigal, g, temp_dir): g for g in genomes_list}
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Gene prediction"):
                faa_files.append(fut.result())
    else:
        faa_files = [g for g in genomes_list if os.path.exists(g)]

    if not faa_files:
        print("No input files found. Exiting.")
        return

    # 2. HMM search
    print("Running hmmsearch in parallel...")
    dbcan_files = parallel_hmmsearch(faa_files, db, temp_dir, num_workers, cpu_per_task=2)

    # 3. Parse HMM output
    print("Parsing hmmsearch output...")
    parsed_files = parallel_parse(dbcan_files, temp_dir, num_workers)

    # 4. Extract features
    print("Extracting features...")
    df_fam, df_sub = extract_features_parallel(parsed_files, subs_dict, num_workers)

    # Optional: save intermediate files
    df_fam.to_csv(os.path.join(temp_dir, 'allfams.tsv'), sep='\t', index=False)
    df_sub.to_csv(os.path.join(temp_dir, 'allsubs.tsv'), sep='\t', index=False)

    # 5. Compute pairwise RPS
    print("Computing pairwise RPS (parallel)...")
    compute_rps_parallel(df_sub, num_workers, output_file)

    # Cleanup (optional – uncomment if you want to delete temp)
    # shutil.rmtree(temp_dir, ignore_errors=True)
    print("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="CaCo parallel: Carbon source competition prediction.")
    parser.add_argument("-m", type=str, default='from_proteins',
                        choices=['download', 'from_proteins', 'from_nucleotides'],
                        help="Mode: download, from_proteins, or from_nucleotides")
    parser.add_argument("-db", type=str, default=None, help="Path to dbCAN HMM database")
    parser.add_argument("-subs", type=str, default=None, help="Path to substrate key JSON")
    parser.add_argument("-tmp", type=str, default='tmp/', help="Temporary directory")
    parser.add_argument("-g1", type=str, help="Genome file 1 (or list file with -gl)")
    parser.add_argument("-g2", type=str, help="Genome file 2")
    parser.add_argument("-gl", type=str, default=None, help="File with list of genome paths")
    parser.add_argument("-o", type=str, default='carboncomp_output.tsv.xz', help="Output file")
    parser.add_argument("-cpus", type=int, default=mp.cpu_count(),
                        help="Number of CPU cores to use (default: all)")
    args = parser.parse_args()

    if args.m == 'download':
        print(download_db())
        sys.exit(0)

    # Determine genome list
    if args.gl:
        with open(args.gl, 'r') as f:
            genomes = [line.strip() for line in f if line.strip()]
    else:
        if args.g1 and args.g2:
            genomes = [args.g1, args.g2]
        else:
            print("Error: Provide either -gl or both -g1 and -g2")
            sys.exit(1)

    # Resolve database and substrate files
    db = args.db if args.db else get_file_path('data/dbcan.hmm')
    subs = args.subs if args.subs else get_file_path('data/substrate_key.json')
    with open(subs, 'r') as f:
        subs_dict = json.load(f)

    main(mode=args.m,
         genomes_list=genomes,
         temp_dir=args.tmp,
         output_file=args.o,
         db=db,
         subs_dict=subs_dict,
         num_workers=args.cpus)
