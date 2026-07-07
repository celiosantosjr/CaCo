# CaCo

[![DBcan Database](https://img.shields.io/badge/DBcan-Database-blue)](https://bcb.unl.edu/dbCAN2/)
[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub issues](https://img.shields.io/github/issues/celiosantosjr/CaCo)](https://github.com/celiosantosjr/CaCo/issues)
[![CI](https://github.com/celiosantosjr/CaCo/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/celiosantosjr/CaCo/actions)

CaCo: A program for predicting carbon source competition and ecological type of interaction from genomes.

To see the recent changes and updates, please follow our [What's New Documentation](whats_new.md).

## Installation

First download the repository:

```
git clone https://github.com/celiosantosjr/CaCo
cd CaCo/
```

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

A. Table of the DBCAN protein families:

    - genome: genome used in the comparison
    - families: list of DBCAN protein families found in the genome

Example:

| genome | families |
| :---: | :---: |
| TARA_SAMN05326651_METAG_RED00102 | AA3, PL12, GH102, CE11, GT2, GT4, GH73, GT19, AA7, GH109, GT5, GT41, CE4, GT9 |
| TARA_SAMEA2620259_METAG_HFKLHEHB | GH18, GT28, CE1, GT30, CE11, GT17, GH1, GT2, GT4, GT9, GH97, GH13, GT5, CE4, GH171, GH2 |


B. Table of the substrates resources:

    - genome: genome used in the comparison
    - substrates: list of substrates found to be possible to be consumed by the genome

Example:

| genome | substrates |
| :---: | :---: |
| TARA_SAMN05326651_METAG_RED00102 | cellooligosaccharide, host glycan, exo-polysaccharide, chitin, chitooligosaccharide, peptidoglycan, cellulose, xylan, glucooligosaccharide, lignin | 
| TARA_SAMEA2620259_METAG_HFKLHEHB | alkaloid, alpha-glucan, alpha-mannan, arabinan, beta-fucosides, beta-galactan, beta-glucan, beta-glucuronan, beta-mannan, chitin, chitosan, exo-polysaccharide, glucosylglycerate, glycogen, host glycan, human milk polysaccharide, pectin, peptidoglycan, polyphenol, starch, sucrose, trehalose, uric acid, xylan | 

C. The table with competition scores contains the following columns:

    - genome1: genome used in the comparison
    - genome2: genome used in the comparsion
    - set1: number of unique carbon substrates detected in genome 1
    - set2: number of unique carbon substrates detected in genome 2
    - intersection: number of carbon substrates overlapping between genome 1 and 2
    - competition: the Jaccard distance between sets 1 and 2
    - relcomp (relative competition): the Jaccard distance between sets 1 and 2 normalized by the maximum Jaccard distance possible
    - prob: Permutation probability calculated with 1,000 simulations of randomly get such intersection, thus competition
    - RPS: Resource Partitioning Score (1 - 2 * competition)
    - relRPS: Resource Partitioning Score (1 - 2 * relative competition)

Example:

| genome1	| genome2	| set1	| set2	| intersection	| competition	| relcomp	| prob	| RPS | relRPS |
| :---: | :---: | :---: | :---: | :---: | :---: |  :---: |  :---: | :---: | :---: |
| TARA_SAMEA2620259_METAG_HFKLHEHB | TARA_SAMN05326651_METAG_RED00102 | 24 | 10 | 5 | 0.172413793 | 0.413793103 | 0.372 | 0.655172414 | 0.172413793 |

