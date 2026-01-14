# CaCo
CaCo: A program for predicting carbon source competition and ecological type of interaction from genomes.

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
| TARA_SAMEA2620259_METAG_HFKLHEHB | GH18, GT28, CE1, GT30, CE11, GT17, GH1, GT2, GT4, GT9, GH97, GH13, GT5, CE4, GH171 |


B. Table of the substrates resources:

    - genome: genome used in the comparison
    - substrates: list of substrates found to be possible to be consumed by the genome

Example:

| genome | substrates |
| :---: | :---: |
| TARA_SAMN05326651_METAG_RED00102 | cellooligosaccharide, host glycan, exo-polysaccharide, chitin, chitooligosaccharide, peptidoglycan, cellulose, xylan, glucooligosaccharide, lignin | 
| TARA_SAMEA2620259_METAG_HFKLHEHB | host glycan, exo-polysaccharide, chitin, glycogen, uric acid, beta-galactan, beta-glucuronan, beta-glucan, glucosylglycerate, xylan, beta-mannan, sucrose, human milk polysaccharide, alpha-glucan, alkaloid, trehalose, polyphenol, peptidoglycan, starch, beta-fucosides | 

C. The table with competition scores contains the following columns:

    - genome1: genome used in the comparison
    - genome2: genome used in the comparsion
    - set1: number of unique carbon substrates detected in genome 1
    - set2: number of unique carbon substrates detected in genome 2
    - intersection: number of carbon substrates overlapping between genome 1 and 2
    - competition: the Jaccard distance between sets 1 and 2 normalized by the maximum Jaccard distance possible
    - prob: Permutation probability calculated with 1,000 simulations of randomly get such intersection, thus competition
    - RPS: Resource Partitioning Score (1 - 2 * competition)

Example:

| genome1	| genome2	| set1	| set2	| intersection	| competition	| prob	| RPS |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| TARA_SAMN05326651_METAG_RED00102 | TARA_SAMEA2620259_METAG_HFKLHEHB | 10 | 20 | 5 | 0.4 | 0.212 | 0.19999999999999996 |
