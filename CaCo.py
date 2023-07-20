import os
import json
import lzma
import math
import argparse
import pandas as pd

from glob import glob
from tqdm import tqdm
from subprocess import run
from itertools import combinations


def get_file_path(infile):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_name = infile  # Replace with your file name
    file_path = os.path.join(script_dir, file_name)
    return file_path


def download_db():
    ofolder = get_file_path('data/')
    os.makedirs(ofolder, exist_ok=True)
    run(['wget',
         '-O', f'{ofolder}/dbcan.hmm',
         'https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V11.txt'])


def dbcansearch(infile, ofolder, db):
    ofile = infile.split("/")[-1].replace(".faa", ".dbcan")
    ofile = f'{ofolder}/{ofile}'
    if os.path.exists(f'{ofile}.parsed'):
        return ofile
        if os.path.exists(ofile):
            return ofile
    else:
        run(['hmmsearch',
             '--cpu', '3',
             '-o', '/dev/null',
             '--domtblout', ofile,
             db,
             infile])
        return ofile
         

def process_hmmsearch_output(infile, ofolder):
    output_file = f'{infile}.parsed'
    if os.path.exists(output_file):
        return ''
    # Define minimum coverage and e-value thresholds
    min_score = 25
    min_coverage = 0.35  # Minimum hit coverage
    max_evalue = 1e-15   # Maximum hit e-value
    # Open the output file from hmmsearch
    with open(infile, 'r') as f:
        lines = f.readlines()
    # Initialize variables to store the best hits
    best_hits = {}
    current_query = None
    current_best_hit = None
    # Iterate over the lines in the output file
    for line in lines:
        # Check if the line contains a query domain
        if line.startswith('#'):
            # Skip the header lines
            continue
        # Extract the fields from the tab-separated line
        fields = line.strip().split()
        # Check if this line corresponds to a new query domain
        query_name = fields[0]
        if query_name != current_query:
            # Store the best hit for the previous query domain, if available
            if current_query is not None and current_best_hit is not None:
                hit_name, hit_evalue, hit_score, hit_start, hit_end, hit_coverage = current_best_hit
                if (hit_coverage >= min_coverage) and (hit_evalue <= max_evalue) and (hit_score >= min_score):
                    best_hits[current_query] = current_best_hit
            # Reset the current query and best hit
            current_query = query_name
            current_best_hit = None
        # Extract hit information from the fields
        hit_name = fields[3].split('.')[0]
        hit_evalue = float(fields[6])
        hit_score = float(fields[7])
        hit_start = int(fields[15])
        hit_end = int(fields[16])
        hit_length = int(fields[17])
        query_length = int(fields[2])
        hit_coverage = (hit_end - hit_start + 1) / query_length
        hit_coverage = min(1.0, hit_coverage)  # fix the coverage when the hit spans across multiple non-contiguous regions of the query sequence.
        # Update the current best hit if needed
        if current_best_hit is None or hit_score > current_best_hit[2] or \
           (hit_score == current_best_hit[2] and hit_evalue < current_best_hit[1]) or \
           (hit_score == current_best_hit[2] and hit_evalue == current_best_hit[1] and hit_coverage > current_best_hit[5]):
            current_best_hit = (hit_name,
                                hit_evalue,
                                hit_score,
                                hit_start,
                                hit_end,
                                hit_coverage)
    # Store the best hit for the last query domain, if available
    if current_query is not None and current_best_hit is not None:
        hit_name, hit_evalue, hit_score, hit_start, hit_end, hit_coverage = current_best_hit
        if hit_coverage >= min_coverage and hit_evalue <= max_evalue:
            best_hits[current_query] = current_best_hit
    # Create a pandas DataFrame from the best hits dictionary
    df = pd.DataFrame.from_dict(best_hits,
                                orient='index',
                                columns=['Hit_Name', 'Hit_Evalue',
                                         'Hit_Score', 'Hit_Start',
                                         'Hit_End', 'Hit_Coverage'])
    # Write the DataFrame to a TSV file
    df.to_csv(output_file,
              sep='\t',
              index_label='Query_Domain')


def EIT(infile, ofile):
    if not ofile.endswith('.xz'): ofile += '.xz'
    infile = pd.read_table(f'{infile}')
    infile = infile.dropna()
    if len(infile) < 2:
        return 'Only genomes without assigned substrates'
    with lzma.open(f'{ofile}', 'wt') as ofile:
        ofile.write('genome1\tgenome2\tset1\tset2\tintersection\tcompetition\trelcomp\tprob\tEIT\trelEIT\n')
        for i, j in combinations(infile.genome, 2):
            try:
                x = calcomp(infile, i, j)
            except:
                print(f'ERR :: Failed the pair {i} and {j}')
            else:
                x = [str(y) for y in x]
                x = '\t'.join(x)
                ofile.write(f'{x}\n')
    
    
def extract_feat(infolder, subs):
    with open(f'allfams.tsv', 'a+') as handle:
        handle.write('genome\tfamilies\n')
        with open(f'allsubs.tsv', 'a+') as handle1:
            handle1.write('genome\tsubstrates\n')
            for infile in tqdm(glob(f'{infolder}/*.dbcan.parsed')):
                namea = infile.split('/')[-1].replace('.dbcan.parsed', '')
                a = pd.read_table(infile)
                a = set(a.Hit_Name)
                a = {x.split('_')[0] for x in a}
                subsa = set(', '.join([subs.get(y, '') for y in a if y in subs]).split(', '))
                handle.write(f"{namea}\t{', '.join(a)}\n")
                handle1.write(f"{namea}\t{', '.join(subsa)}\n")


def calculate_overlap_probability(n1, n2, k, num_trials=1000):
    import random
    random.seed(42)
    M = 60  # max substrates possible
    overlap_count = 0
    for _ in range(num_trials):
        sample1 = random.sample(range(M), min(n1, M))
        sample2 = random.sample(range(M), min(n2, M))
        overlap = len(set(sample1).intersection(set(sample2)))
        if overlap >= k:
            overlap_count += 1
    probability = overlap_count / num_trials
    return probability
    
    
def probcomp(m, n, k):
    # initial check
    if (m+n) <= k:
        print('ERR :: Intersection is wrong sized')
        return 'n.a.', 'n.a.', 'n.a.'
    if (m == 0) and (n == 0):
        print('ERR :: Both sets are empty')
        return 0, 0, 1.0
    if (m == 0) or (n == 0):
        print('WARN :: A set is empty')
        return 0, 0, 1.0
    if not ((m >= k) and (n >= k)):
        print('ERR :: Intersection is bigger than ind. set')
        return 'n.a.', 'n.a.'
    # to calculate maxcomp
    ke = min(m, n)
    maxcomp = ke / (m+n-ke)
    comp = k/(m+n-k)
    comp /= maxcomp
    probability = calculate_overlap_probability(m, n, k)
    return k/(m+n-k), comp, probability
    
    
def calcomp(df, x1, x2):
    a = set(df.loc[df.genome==x1][df.columns[-1]].tolist()[0].split(', '))
    b = set(df.loc[df.genome==x2][df.columns[-1]].tolist()[0].split(', '))
    intersect = len(a.intersection(b))
    competition, relcomp, prob = probcomp(len(a), len(b), intersect)
    EIT = 1 - 2*competition
    relEIT = 1 - 2*relcomp
    return (x1, x2, len(a),
            len(b), intersect,
            competition, relcomp, prob,
            EIT, relEIT)    


def predict_genes(infile, outdir):
    ofile = infile.split('/')[-1]
    ofile = '.'.join(ofile.split('.')[:-1])
    ofile += '.faa'
    run(['prodigal',
         '-a', f'{outdir}/{ofile}',
         '-c',
         '-i', infile,
         '-m',
         '-o', '/dev/null',
         '-p', 'single',
         '-q'])
         
         
def main(mode=None, ofile=None, temp=None, genomes=None, g1=None, g2=None, db=None, subs=None):
    subs = json.load(open(subs))
    os.makedirs(temp, exist_ok=True)
    if genomes:
            genomes = [row.strip() for row in open(genomes, 'r')]
            
    def process_input_files(inputs):
        for infile in inputs:
            infile = dbcansearch(infile, temp, db)
            process_hmmsearch_output(infile, temp)
            run(['rm', '-rf', infile])
            print(f'\t{infile}...')
    
    if mode == 'from_nucleotides':
        print('Gene prediction')
        inputs = genomes or [g1, g2]
        for g in inputs:
            predict_genes(g, temp)
            print(f'\t{g}...')

        print('Annotation using DBCan')
        process_input_files(glob(f'{temp}/*.faa'))
    else:
        print('Annotation using DBCan')
        inputs = genomes or [g1, g2]
        process_input_files(inputs)
    
    print('Extracting features')
    extract_feat(temp, subs)
    
    print('Calculating EIT for substrates')
    EIT('allsubs.tsv', ofile)

    print('Cleaning')
    run(['rm', '-rf', temp])
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="CaCo: A program for predicting carbon source competition and ecological type of interaction from genomes.")
    parser.add_argument("-m", type=str, default='from_proteins', required=True, help="Mode: download, from_proteins or from_nucleotides (default: from_proteins)")
    parser.add_argument("-db", type=str, default=None, help="Address to the dbCAN database")
    parser.add_argument("-subs", type=str, default=None, help="Address to the substrate database")
    parser.add_argument("-tmp", type=str, default='tmp/', help='Temporary directory address (default: tmp/)')
    parser.add_argument("-g1", type=str, help="Address to genome 1")
    parser.add_argument("-g2", type=str, help="Address to genome 2")
    parser.add_argument("-gl", type=str, default=None, help="List of genome addresses to test pairwise comparison override g1 and g2 arguments")
    parser.add_argument("-o", type=str, default='carboncomp_output.tsv.xz', help='Output file (default: carboncomp_output.tsv.xz)')
       
    args = parser.parse_args()

    mode = args.m
    
    if mode == 'download':
        download_db()
        print('Downloaded DBCAN')
    else:
        if args.gl == None:
            g1 = args.g1
            g2 = args.g2
        else:
            gl = args.gl    
        if args.db == None:
            database=get_file_path('data/dbcan.hmm')
        else:
            database = args.db
        if args.subs == None:
            subs=get_file_path('data/substrate_key.json')
        else:
            subs = args.subs
        output_file = args.o
        temp = args.tmp
        if args.gl:
            main(mode=mode,
                 ofile=output_file,
                 temp=temp,
                 db=database,
                 subs=subs,
                 genomes=gl)
        else:
            main(mode=mode,
                 ofile=output_file,
                 temp=temp,
                 db=database,
                 subs=subs,
                 g1=g1,
                 g2=g2)
    
