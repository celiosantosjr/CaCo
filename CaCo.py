#!/usr/bin/env python3
# CaCo.py - High-throughput Carbon Competition prediction

import os
import sys
import json
import lzma
import argparse
import tempfile
import shutil
import random
import subprocess
import multiprocessing as mp
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm
try:
    import pyrodigal
except ImportError:
    pyrodigal = None


# ---------- Version ------------------------------------------------
__version__ = "0.2.1a"

# ---------- Helper: get script directory (for data files) ----------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def get_data_path(infile):
    """Resolve paths relative to script directory for static data."""
    return os.path.join(SCRIPT_DIR, infile)


# ---------- Download dbCAN HMM database ----------
def download_db():
    ofolder = get_data_path('data/')
    os.makedirs(ofolder, exist_ok=True)
    file_path = os.path.join(ofolder, 'dbcan.hmm')
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


# ---------- Gene prediction: two backends ----------
def predict_genes_prodigal(infile, outdir):
    """Use original Prodigal (subprocess) for exact reproducibility."""
    outfile = os.path.join(outdir, os.path.basename(infile).rsplit('.', 1)[0] + '.faa')
    if os.path.exists(outfile):
        return outfile
    cmd = ['prodigal', '-i', infile, '-a', outfile, '-p', 'meta', '-q']
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return outfile

def predict_genes_pyrodigal(infile, outdir):
    """Fast Python binding – may produce slightly different calls."""
    if pyrodigal is None:
        raise RuntimeError("Pyrodigal is required for --use-pyrodigal. Install it or omit that option.")
    outfile = os.path.join(outdir, os.path.basename(infile).rsplit('.', 1)[0] + '.faa')
    if os.path.exists(outfile):
        return outfile
    with open(infile, 'rb') as f:
        seq = f.read()
    genes = pyrodigal.GeneFinder(meta=False)
    genes.train(seq)
    with open(outfile, 'w') as out:
        for idx, pred in enumerate(genes.find_genes(seq), 1):
            prot = pred.translate()
            out.write(f">{os.path.basename(infile)}_gene_{idx}\n{prot}\n")
    return outfile


def predict_genes_parallel(genomes_list, outdir, num_workers, use_pyrodigal):
    """Dispatch gene prediction in parallel using the chosen backend."""
    predictor = predict_genes_pyrodigal if use_pyrodigal else predict_genes_prodigal
    if not use_pyrodigal:
        # Check that Prodigal is installed
        if shutil.which('prodigal') is None:
            raise RuntimeError("Prodigal not found in PATH. Install it or use --use-pyrodigal.")
    faa_files = []
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = {ex.submit(predictor, g, outdir): g for g in genomes_list}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Gene prediction"):
            faa_files.append(fut.result())
    return faa_files


# ---------- HMM search with hmmsearch (parallel) ----------
def run_hmmsearch(genome_faa, db, outdir, cpu_per_task):
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


def parallel_hmmsearch(faa_files, db, outdir, total_cpus):
    """
    Run hmmsearch in parallel without oversubscribing.
    We compute: num_jobs = min(total_cpus, len(faa_files))
    cpu_per_task = max(1, total_cpus // num_jobs)
    This ensures total requested threads <= total_cpus.
    """
    num_jobs = min(total_cpus, len(faa_files))
    cpu_per_task = max(1, total_cpus // num_jobs)
    with ProcessPoolExecutor(max_workers=num_jobs) as ex:
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


# ---------- Feature extraction (parallel, memory‑safe) ----------
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


def matrix_to_feature_tables(matrix_file, output_dir, id_column=None):
    """Read a genome-by-annotation table and convert positive values to presence/absence.

    The input may be tab-delimited or comma-separated. The first column is used as
    the genome identifier by default; use id_column to select a named identifier
    column. Every remaining column is treated as an annotation/substrate feature.
    Values greater than zero are converted to 1, while zero, missing, and nonnumeric
    values are converted to 0. The resulting substrate lists are compatible with
    the existing pairwise RPS calculation.
    """
    if not os.path.isfile(matrix_file):
        raise FileNotFoundError(f"Matrix file not found: {matrix_file}")

    try:
        df = pd.read_csv(matrix_file, sep=None, engine='python')
    except Exception as e:
        raise ValueError(f"Could not read matrix as CSV/TSV: {e}") from e

    if df.empty or df.shape[1] < 2:
        raise ValueError("The matrix must contain an identifier column and at least one annotation column.")

    if id_column is None:
        id_column = df.columns[0]
    if id_column not in df.columns:
        raise ValueError(f"Identifier column '{id_column}' was not found. Available columns: {list(df.columns)}")

    if df[id_column].isna().any() or (df[id_column].astype(str).str.strip() == '').any():
        raise ValueError(f"Identifier column '{id_column}' contains missing or empty genome identifiers.")
    if df[id_column].duplicated().any():
        duplicates = df.loc[df[id_column].duplicated(), id_column].astype(str).tolist()
        raise ValueError(f"Genome identifiers must be unique; duplicated values include: {duplicates[:5]}")

    feature_columns = [c for c in df.columns if c != id_column]
    numeric = df[feature_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    binary = (numeric > 0).astype('uint8')

    records = []
    for genome, row in zip(df[id_column].astype(str), binary.itertuples(index=False, name=None)):
        present = [feature for feature, value in zip(feature_columns, row) if value == 1]
        records.append({
            'genome': genome,
            'families': ', '.join(present),
            'substrates': ', '.join(present)
        })

    fam_path = os.path.join(output_dir, 'allfams.tsv')
    sub_path = os.path.join(output_dir, 'allsubs.tsv')
    pd.DataFrame([{'genome': r['genome'], 'families': r['families']} for r in records]).to_csv(
        fam_path, sep='\t', index=False
    )
    pd.DataFrame([{'genome': r['genome'], 'substrates': r['substrates']} for r in records]).to_csv(
        sub_path, sep='\t', index=False
    )
    return fam_path, sub_path


def extract_features_parallel(parsed_files, subs_dict, num_workers, output_dir):
    """
    Extract features in parallel and write results to disk.

    Note: Results are accumulated in memory (as DataFrames) and written once
    at the end. This is fast for up to ~10,000 genomes. For larger datasets,
    a streaming/chunked writer would be required to avoid memory exhaustion.
    """
    fam_path = os.path.join(output_dir, 'allfams.tsv')
    sub_path = os.path.join(output_dir, 'allsubs.tsv')

    if len(parsed_files) > 10000:
        print("Warning: many genomes – consider a chunked writer to avoid memory issues.")
    fam_records = []
    sub_records = []
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = {ex.submit(extract_features_one, pf, subs_dict): pf for pf in parsed_files}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Extracting features"):
            name, fams, subs = fut.result()
            fam_records.append({'genome': name, 'families': ', '.join(fams)})
            sub_records.append({'genome': name, 'substrates': ', '.join(subs)})
    # Write out
    pd.DataFrame(fam_records).to_csv(fam_path, sep='\t', index=False)
    pd.DataFrame(sub_records).to_csv(sub_path, sep='\t', index=False)
    return fam_path, sub_path


# ---------- Probability overlap ----------
def probability_overlap(n1, n2, k, trials=1000):
    random.seed(42)
    M = 60
    # Use a hash-based seed that depends on the unordered pair
    # This ensures the same result regardless of order
    seed = hash((min(n1, n2), max(n1, n2), k, M))
    rng = random.Random(seed)
    
    a = min(n1, M)
    b = min(n2, M)
    
    # Ensure a <= b for consistent sampling order
    if a > b:
        a, b = b, a
    
    cnt = 0
    for _ in range(trials):
        s1 = set(rng.sample(range(M), a))
        s2 = set(rng.sample(range(M), b))
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
def compute_rps_chunk(pair_chunk, df_sub, temp_dir, chunk_id):
    temp_out = os.path.join(temp_dir, f"rps_{os.getpid()}_{chunk_id}.tsv")
    with open(temp_out, 'w') as f:
        f.write("genome1\tgenome2\tset1\tset2\tintersection\tcompetition\trelcomp\tprob\tRPS\trelRPS\n")
        for g1, g2 in pair_chunk:
            # Ensure deterministic order: smaller genome name first
            if g1 > g2:
                g1, g2 = g2, g1
            row1 = df_sub[df_sub.genome == g1]
            row2 = df_sub[df_sub.genome == g2]
            if row1.empty or row2.empty:
                continue
            substrate_text1 = row1.iloc[0]['substrates']
            substrate_text2 = row2.iloc[0]['substrates']
            set1 = set() if pd.isna(substrate_text1) or not str(substrate_text1) else set(str(substrate_text1).split(', '))
            set2 = set() if pd.isna(substrate_text2) or not str(substrate_text2) else set(str(substrate_text2).split(', '))
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
            futures = [ex.submit(compute_rps_chunk, ch, df_sub, temp_dir, chunk_id)
                       for chunk_id, ch in enumerate(chunks)]
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
def main(mode, genomes_list, temp_dir, output_file, db, subs_dict, num_workers, use_pyrodigal,
         matrix_file=None, matrix_id_column=None):
    # Determine output directory: all final outputs go to current working directory
    output_dir = os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)

    # Matrix mode: annotations are treated as substrate features after binarization.
    if mode == 'from_matrix':
        print("Reading annotation matrix and converting positive values to binary presence...")
        fam_path, sub_path = matrix_to_feature_tables(matrix_file, output_dir, matrix_id_column)
    # 1. Gene prediction (if nucleotides)
    elif mode == 'from_nucleotides':
        print("Gene prediction...")
        faa_files = predict_genes_parallel(genomes_list, temp_dir, num_workers, use_pyrodigal)
    else:
        faa_files = [g for g in genomes_list if os.path.exists(g)]

    if mode != 'from_matrix' and not faa_files:
        print("No input files found. Exiting.")
        return

    if mode == 'from_matrix':
        print(f"Reading substrate table: {sub_path}")
    else:
        # 2. HMM search (with dynamic CPU allocation to avoid oversubscription)
        print("Running hmmsearch in parallel...")
        dbcan_files = parallel_hmmsearch(faa_files, db, temp_dir, num_workers)

        # 3. Parse HMM output
        print("Parsing hmmsearch output...")
        parsed_files = parallel_parse(dbcan_files, temp_dir, num_workers)

        # 4. Extract features (write to output_dir)
        print("Extracting features...")
        fam_path, sub_path = extract_features_parallel(parsed_files, subs_dict, num_workers, output_dir)

    # 5. Read the substrate table for RPS computation
    df_sub = pd.read_table(sub_path)

    # 6. Compute pairwise RPS (output to output_dir)
    out_path = os.path.join(output_dir, output_file) if not os.path.isabs(output_file) else output_file
    print("Computing pairwise RPS (parallel)...")
    compute_rps_parallel(df_sub, num_workers, out_path)

    print(f"Done. Output files: {fam_path}, {sub_path}, {out_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="CaCo parallel: Carbon source competition prediction.")
    parser.add_argument("-v", "--version", action="version", version=f"CaCo {__version__}")
    parser.add_argument("-m", type=str, default='from_proteins',
                        choices=['download', 'from_proteins', 'from_nucleotides', 'from_matrix'],
                        help="Mode: download, from_proteins, from_nucleotides, or from_matrix")
    parser.add_argument("-db", type=str, default=None, help="Path to dbCAN HMM database")
    parser.add_argument("-subs", type=str, default=None, help="Path to substrate key JSON")
    parser.add_argument("-tmp", type=str, default='tmp/', help="Temporary directory")
    parser.add_argument("-g1", type=str, help="Genome file 1 (or list file with -gl)")
    parser.add_argument("-g2", type=str, help="Genome file 2")
    parser.add_argument("-gl", type=str, default=None, help="File with list of genome paths")
    parser.add_argument("-matrix", type=str, default=None,
                        help="Genome-by-annotation CSV/TSV table for from_matrix mode")
    parser.add_argument("--matrix-id-column", type=str, default=None,
                        help="Column containing genome identifiers; defaults to the first column")
    parser.add_argument("-o", type=str, default='carboncomp_output.tsv.xz', help="Output file name (placed in current directory)")
    parser.add_argument("-cpus", type=int, default=mp.cpu_count(),
                        help="Number of CPU cores to use (default: all)")
    parser.add_argument("--use-pyrodigal", action="store_true",
                        help="Use Pyrodigal for gene prediction (faster, but slightly different from the original PNAS code). If not set, uses Prodigal (subprocess) for exact reproducibility.")
    args = parser.parse_args()

    if args.m == 'download':
        print(download_db())
        sys.exit(0)

    if args.m == 'from_matrix':
        if not args.matrix:
            print("Error: -matrix is required when using -m from_matrix")
            sys.exit(1)
        genomes = []
    # Determine genome list for FASTA/protein modes
    elif args.gl:
        with open(args.gl, 'r') as f:
            genomes = [line.strip() for line in f if line.strip()]
    else:
        if args.g1 and args.g2:
            genomes = [args.g1, args.g2]
        else:
            print("Error: Provide either -gl or both -g1 and -g2")
            sys.exit(1)

    # Resolve database and substrate files. Matrix mode does not need either file.
    if args.m == 'from_matrix':
        db = None
        subs_dict = {}
    else:
        db = args.db if args.db else get_data_path('data/dbcan.hmm')
        subs = args.subs if args.subs else get_data_path('data/substrate_key.json')
        with open(subs, 'r') as f:
            subs_dict = json.load(f)

    if args.m == 'from_matrix':
        print("Using matrix mode: all positive annotation values will be treated as present.")
    elif args.use_pyrodigal:
        print("Using Pyrodigal for gene prediction (fast). Note: results may differ slightly from the original paper.")
    else:
        print("Using Prodigal (subprocess) for exact reproducibility with the PNAS paper.")

    main(mode=args.m,
         genomes_list=genomes,
         temp_dir=args.tmp,
         output_file=args.o,
         db=db,
         subs_dict=subs_dict,
         num_workers=args.cpus,
         use_pyrodigal=args.use_pyrodigal,
         matrix_file=args.matrix,
         matrix_id_column=args.matrix_id_column)
