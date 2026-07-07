import os
import sys
import json
import lzma
import math
import argparse
import tempfile
import shutil
from glob import glob
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import pandas as pd
from tqdm import tqdm
import pyrodigal  # Fast Prodigal replacement


# ---------- Helper functions ----------
def get_file_path(infile):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, infile)


def download_db():
    """Download dbCAN HMM database (unchanged)."""
    ofolder = get_file_path('data/')
    os.makedirs(ofolder, exist_ok=True)
    file_path = f'{ofolder}/dbcan.hmm'
    from subprocess import run
    result = run(['wget',
                  '-O', file_path,
                  'https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-V11.txt'],
                 capture_output=True, text=True)
    if result.returncode != 0:
        return f"Download failed: {result.stderr}"
    # Validate file header and footer (simplified)
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


# ---------- Gene prediction with Pyrodigal ----------
def predict_genes_pyrodigal(infile, outdir, threads=1):
    """Run Pyrodigal on a single genome; returns output file path."""
    outfile = os.path.join(outdir, os.path.basename(infile).rsplit('.', 1)[0] + '.faa')
    # Use the Pyrodigal ORF finder
    with open(infile, 'rb') as f:
        sequence = f.read()
    # Train and predict genes (single mode, meta mode can be set via flag)
    # We keep -p single as in original
    genes = pyrodigal.GeneFinder(meta=False)  # single genome mode
    genes.train(sequence)
    # Write protein sequences in FASTA format
    with open(outfile, 'w') as out:
        for idx, pred in enumerate(genes.find_genes(sequence), 1):
            # pred is a Gene object with attributes: start, end, strand, sequence
            # We need to get the protein sequence; pyrodigal provides .translate()
            protein = pred.translate()
            header = f">{os.path.basename(infile)}_gene_{idx}"
            out.write(f"{header}\n{protein}\n")
    return outfile


# ---------- Parallel HMM search with hmmsearch (using multiprocessing) ----------
def run_hmmsearch(genome_faa, db, outdir, cpu_per_task=2):
    """Run hmmsearch on one FASTA file and return path to domtblout."""
    outfile = os.path.join(outdir, os.path.basename(genome_faa).replace('.faa', '.dbcan'))
    if os.path.exists(outfile):
        return outfile  # already done
    from subprocess import run
    cmd = ['hmmsearch',
           '--cpu', str(cpu_per_task),
           '-o', '/dev/null',
           '--domtblout', outfile,
           db, genome_faa]
    run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return outfile


def parallel_hmmsearch(faa_files, db, outdir, num_workers, cpu_per_task=2):
    """Launch hmmsearch for all FASTA files in parallel."""
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(run_hmmsearch, faa, db, outdir, cpu_per_task): faa
                   for faa in faa_files}
        results = []
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc="HMM search"):
            results.append(future.result())
    return results


# ---------- Parse hmmsearch output (parallel) ----------
def parse_hmmsearch_output(infile, outdir):
    """Process a single .dbcan file and write parsed .parsed file."""
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
            query_name = fields[0]
            hit_name = fields[3].split('.')[0]
            hit_evalue = float(fields[6])
            hit_score = float(fields[7])
            hit_start = int(fields[15])
            hit_end = int(fields[16])
            query_len = int(fields[2])
            coverage = min(1.0, (hit_end - hit_start + 1) / query_len)

            if query_name != current_query:
                if current_query is not None and current_best is not None:
                    if (current_best[5] >= min_coverage and
                        current_best[1] <= max_evalue and
                        current_best[2] >= min_score):
                        best_hits[current_query] = current_best
                current_query = query_name
                current_best = None

            if (current_best is None or hit_score > current_best[2] or
                (hit_score == current_best[2] and hit_evalue < current_best[1]) or
                (hit_score == current_best[2] and hit_evalue == current_best[1] and coverage > current_best[5])):
                current_best = (hit_name, hit_evalue, hit_score, hit_start, hit_end, coverage)

    # Last query
    if current_query is not None and current_best is not None:
        if (current_best[5] >= min_coverage and
            current_best[1] <= max_evalue and
            current_best[2] >= min_score):
            best_hits[current_query] = current_best

    # Write to parsed file
    df = pd.DataFrame.from_dict(best_hits, orient='index',
                                columns=['Hit_Name','Hit_Evalue','Hit_Score',
                                         'Hit_Start','Hit_End','Hit_Coverage'])
    df.to_csv(parsed_file, sep='\t', index_label='Query_Domain')
    return parsed_file


def parallel_parse(dbcan_files, outdir, num_workers):
    """Parse all .dbcan files in parallel."""
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(parse_hmmsearch_output, f, outdir): f
                   for f in dbcan_files}
        results = []
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc="Parsing HMM output"):
            results.append(future.result())
    return results


# ---------- Feature extraction (parallel) ----------
def extract_features_for_one(parsed_file, subs_dict):
    """Extract families and substrates from one parsed file."""
    name = os.path.basename(parsed_file).replace('.dbcan.parsed', '')
    df = pd.read_table(parsed_file)
    families = set(df.Hit_Name)
    families = {x.split('_')[0] for x in families}   # family IDs
    sorted_fam = sorted(families)
    # substrates
    subs_set = set()
    for fam in families:
        if fam in subs_dict:
            subs_set.update(subs_dict[fam].split(', '))
    sorted_subs = sorted(subs_set)
    return (name, sorted_fam, sorted_subs)


def extract_features_parallel(parsed_files, subs_dict, num_workers):
    """Extract features from all parsed files in parallel, returning two DataFrames."""
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(extract_features_for_one, pf, subs_dict): pf
                   for pf in parsed_files}
        results = []
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc="Extracting features"):
            results.append(future.result())
    # Build DataFrames
    fam_records = []
    sub_records = []
    for name, fams, subs in results:
        fam_records.append({'genome': name, 'families': ', '.join(fams)})
        sub_records.append({'genome': name, 'substrates': ', '.join(subs)})
    df_fam = pd.DataFrame(fam_records)
    df_sub = pd.DataFrame(sub_records)
    return df_fam, df_sub


# ---------- RPS calculation (parallel pairwise) ----------
def probability_overlap(n1, n2, k, trials=1000):
    """Empirical probability of overlap >= k under random sampling."""
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
    """Compute competition and RPS for two sets."""
    if (m+n) <= k or (m==0 and n==0):
        return (0, 0, 1.0)
    if m==0 or n==0:
        return (0, 0, 1.0)
    if not (m >= k and n >= k):
        return ('n.a.', 'n.a.', 'n.a.')
    maxcomp = min(m,n) / (m+n-min(m,n))
    comp = k/(m+n-k)
    relcomp = comp / maxcomp if maxcomp != 0 else 0
    prob = probability_overlap(m, n, k)
    return comp, relcomp, prob


def compute_rps_pair(args):
    """Worker function for one pair (genome1, genome2, df_sub)."""
    idx, (g1, g2, df_sub) = args
    row1 = df_sub[df_sub.genome == g1]
    row2 = df_sub[df_sub.genome == g2]
    if row1.empty or row2.empty:
        return None
    set1 = set(row1.iloc[0]['substrates'].split(', '))
    set2 = set(row2.iloc[0]['substrates'].split(', '))
    intersect = len(set1 & set2)
    m, n = len(set1), len(set2)
    competition, relcomp, prob = probcomp(m, n, intersect)
    if isinstance(competition, str):
        return None
    RPS = 1 - 2 * competition
    relRPS = 1 - 2 * relcomp
    return (g1, g2, m, n, intersect, competition, relcomp, prob, RPS, relRPS)


def compute_rps_parallel(df_sub, num_workers, output_file):
    """Compute RPS for all genome pairs in parallel, writing to output."""
    genomes = df_sub['genome'].tolist()
    pairs = list(combinations(genomes, 2))
    if not pairs:
        return

    # Prepare arguments for workers: we pass each pair and the whole df
    # But df is large; better to pass indices? To avoid pickling large df repeatedly,
    # we can create a list of args with (g1, g2, df_sub) but that copies df each time.
    # We'll use a shared memory approach? Simpler: pass the entire df to each worker via
    # initializer? We can use multiprocessing.Pool with initializer to share df.
    # But ProcessPoolExecutor doesn't support initializer easily. Use multiprocessing.Pool.
    # We'll create a Pool with initializer to set global df_sub.
    # However, to simplify, we'll just pass df_sub as a global variable.
    # We'll define a global variable and use Pool.

    # We'll use multiprocessing.Pool with initializer.
    manager = mp.Manager()
    # Create a shared dict? Actually we can just rely on copy-on-write? For large df, copy may be heavy.
    # Better: use a global variable in the worker function (defined at module level).
    # We'll set a global variable before calling pool.

    # To avoid pickling, we can define a worker that reads df from a global.
    # We'll set a global variable _DF_SUB in the main process, and then use Pool with initializer to set it in each worker.
    # But multiprocessing.Pool initializer cannot easily set global in child; we can use a function that sets it.
    # Simpler: just pass the df as argument; it will be pickled once per worker? Actually it will be pickled for each call.
    # That's inefficient. Better: use a shared memory via Manager.dict? Not needed.

    # Given the number of pairs may be large, we can chunk them and write per worker to temp files.
    # We'll create a Pool with a custom initializer to load df_sub from a global variable.
    # We'll define a global variable DF_SUB and set it before pool creation.

    # We'll set a global variable in the main module.
    global DF_SUB
    DF_SUB = df_sub  # noqa

    def worker(pair):
        g1, g2 = pair
        row1 = DF_SUB[DF_SUB.genome == g1]
        row2 = DF_SUB[DF_SUB.genome == g2]
        if row1.empty or row2.empty:
            return None
        set1 = set(row1.iloc[0]['substrates'].split(', '))
        set2 = set(row2.iloc[0]['substrates'].split(', '))
        intersect = len(set1 & set2)
        m, n = len(set1), len(set2)
        competition, relcomp, prob = probcomp(m, n, intersect)
        if isinstance(competition, str):
            return None
        RPS = 1 - 2 * competition
        relRPS = 1 - 2 * relcomp
        return (g1, g2, m, n, intersect, competition, relcomp, prob, RPS, relRPS)

    # Chunk pairs to write to separate temp files to avoid locking
    chunk_size = max(1, len(pairs) // num_workers)
    chunks = [pairs[i:i+chunk_size] for i in range(0, len(pairs), chunk_size)]

    # We'll use a pool to process each chunk in a worker, and that worker writes to its own temp file.
    def process_chunk(chunk, temp_dir):
        """Process a list of pairs and write results to a temp file."""
        temp_out = os.path.join(temp_dir, f"rps_{os.getpid()}.tsv")
        with open(temp_out, 'w') as f:
            # Write header
            f.write("genome1\tgenome2\tset1\tset2\tintersection\tcompetition\trelcomp\tprob\tRPS\trelRPS\n")
            for g1, g2 in chunk:
                row1 = DF_SUB[DF_SUB.genome == g1]
                row2 = DF_SUB[DF_SUB.genome == g2]
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

    # Create a temporary directory for chunk results
    temp_dir = tempfile.mkdtemp(prefix="rps_chunks_")
    try:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(process_chunk, chunk, temp_dir) for chunk in chunks]
            temp_files = []
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="Computing RPS pairs"):
                temp_files.append(future.result())

        # Merge all temp files into final compressed output
        with lzma.open(output_file, 'wt') as out_f:
            header_written = False
            for tf in temp_files:
                with open(tf, 'r') as inf:
                    if not header_written:
                        # Write header once
                        header = inf.readline()
                        out_f.write(header)
                        header_written = True
                    else:
                        # Skip header in subsequent files
                        inf.readline()
                    # Write data
                    for line in inf:
                        out_f.write(line)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


# ---------- Main workflow ----------
def main(mode, genomes_list, temp_dir, output_file, db, subs_dict, num_workers):
    os.makedirs(temp_dir, exist_ok=True)

    # Step 1: Gene prediction if needed
    if mode == 'from_nucleotides':
        print("Gene prediction using Pyrodigal (parallel)...")
        faa_files = []
        # We can parallelize gene prediction too
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(predict_genes_pyrodigal, g, temp_dir): g
                       for g in genomes_list}
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="Gene prediction"):
                faa_files.append(future.result())
    else:
        # mode == 'from_proteins' – assume genomes_list are already .faa files
        faa_files = [g for g in genomes_list if os.path.exists(g)]

    if not faa_files:
        print("No input files found. Exiting.")
        return

    # Step 2: HMM search (parallel)
    print("Running hmmsearch in parallel...")
    dbcan_files = parallel_hmmsearch(faa_files, db, temp_dir, num_workers,
                                     cpu_per_task=2)  # each hmmsearch uses 2 threads

    # Step 3: Parse hmmsearch output (parallel)
    print("Parsing hmmsearch output in parallel...")
    parsed_files = parallel_parse(dbcan_files, temp_dir, num_workers)

    # Step 4: Extract features (parallel)
    print("Extracting features...")
    df_fam, df_sub = extract_features_parallel(parsed_files, subs_dict, num_workers)

    # Write feature files (optional)
    df_fam.to_csv(os.path.join(temp_dir, 'allfams.tsv'), sep='\t', index=False)
    df_sub.to_csv(os.path.join(temp_dir, 'allsubs.tsv'), sep='\t', index=False)

    # Step 5: Compute RPS pairwise (parallel with chunking)
    print("Computing RPS for all pairs (parallel)...")
    compute_rps_parallel(df_sub, num_workers, output_file)

    # Cleanup temp if not needed (optional)
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

    # Determine list of genomes
    if args.gl:
        with open(args.gl, 'r') as f:
            genomes = [line.strip() for line in f if line.strip()]
    else:
        if args.g1 and args.g2:
            genomes = [args.g1, args.g2]
        else:
            print("Error: Provide either -gl or both -g1 and -g2")
            sys.exit(1)

    # Set paths
    if args.db is None:
        db = get_file_path('data/dbcan.hmm')
    else:
        db = args.db
    if args.subs is None:
        subs = get_file_path('data/substrate_key.json')
    else:
        subs = args.subs

    with open(subs, 'r') as f:
        subs_dict = json.load(f)

    main(mode=args.m,
         genomes_list=genomes,
         temp_dir=args.tmp,
         output_file=args.o,
         db=db,
         subs_dict=subs_dict,
         num_workers=args.cpus)
