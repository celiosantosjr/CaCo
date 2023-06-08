# CaCo
CaCo: A program for predicting carbon source competition and ecological type of interaction from genomes.

## Installation

CaCo is easily installed using conda:

```
conda env create -f CaCo_env.yml
conda activate CaCo
python3 CaCo.py -h
```

Then, you need to download the database from DBCan resource, in case you still do not have it locally.
CaCo provides a simple way to do it:

```
python3 CaCo.py -m download
```

## Usage

CaCo is a program that allows the user to carry on comparisons across groups of genomes or pairs of genomes
easily. To compare pairs you can simply add the address to each one of the genomes and let CaCo work, otherwise
you can add a list of genome addresses and CaCo will compare them all pairwise. CaCo accepts nucleotides and
gene prediction outputs as predicted proteins.

**Do not input genes, only contigs or predicted proteins**

Examples of usage include:

1. Genome groups in the file:

```
python3 CaCo.py -m from_nucleotides -gl example/gltest.txt 
```

2. Genome pairs as contigs:

```
python3 CaCo.py -m from_nucleotides -g1 example/TARA_SAMEA2620259_METAG_HFKLHEHB.fa -g2 example/TARA_SAMN05326651_METAG_RED00102.fa
```

3. Genome pairs as proteins:

```
python3 CaCo.py -m from_proteins -g1 example/TARA_SAMEA2620259_METAG_HFKLHEHB.faa -g2 example/TARA_SAMN05326651_METAG_RED00102.faa
```

To carefully check the options, type:

`python3 CaCo.py -h`

## Output examples

CaCo produces 3 files:

A.
B. 
C. The table with competition scores
