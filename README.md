# CaCo

[![DBcan Database](https://img.shields.io/badge/DBcan-Database-blue)](https://bcb.unl.edu/dbCAN2/)
[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub issues](https://img.shields.io/badge/Issues-Open-red)](https://github.com/celiosantosjr/CaCo/issues)
[![CI](https://github.com/celiosantosjr/CaCo/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/celiosantosjr/CaCo/actions)

**CaCo** is a high-throughput program for predicting carbon-source competition and ecological interaction types from microbial genomes or precomputed genome-by-annotation tables.

> **New in this version:** CaCo supports three analysis inputs: nucleotide FASTA files, predicted protein FASTA files, and genome-by-annotation CSV/TSV tables. The workflow is fully parallelized, uses dynamic CPU allocation, and provides dual-backend gene prediction through Prodigal or Pyrodigal.
>
> See [What's New](whats_new.md) for performance changes and other updates.

---

## Important note on reproducibility

The version of CaCo used in the PNAS paper, “Resource availability structures microbial competition through genomic niche partitioning,” by Santos-Júnior et al., used **Prodigal** through a subprocess for gene prediction.

For FASTA workflows, the default behavior remains Prodigal in parallel. This preserves compatibility with the original gene-calling workflow. Users may add `--use-pyrodigal` to use the faster Pyrodigal Python binding; this can produce slightly different gene calls in some cases and therefore minor downstream score differences.

> For exact reproduction of the PNAS results, run FASTA analyses without `--use-pyrodigal`, or use the [archived commit](https://github.com/celiosantosjr/CaCo/blob/db50a80e974b4daedb217cf971a68c6079fc33e7/CaCo.py) from January 14th.

Matrix mode does not perform gene prediction or dbCAN searches. It assumes that annotation features have already been computed and uses them directly as binary resource features after conversion.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/celiosantosjr/CaCo
cd CaCo/
```

Create and activate the Conda environment:

```bash
conda env create -f CaCo_env.yml
conda activate CaCo
```

The environment installs the dependencies required for FASTA workflows, including **Prodigal** and **Pyrodigal**. Matrix mode requires the Python dependencies used by CaCo, including pandas, but does not require Prodigal, Pyrodigal, the dbCAN HMM database, or the substrate-key JSON file.

If Prodigal needs to be installed separately, use:

```bash
conda install -c bioconda prodigal
```

For nucleotide and protein FASTA workflows, download the dbCAN HMM database with:

```bash
python3 CaCo.py -m download
```

Check the command-line interface with:

```bash
python3 CaCo.py -h
```

---

## Input modes

CaCo supports the following modes:

| Mode | Input | Main processing |
| :--- | :--- | :--- |
| `from_nucleotides` | Complete nucleotide genome or contig FASTA files | Gene prediction, dbCAN HMM search, substrate mapping, and RPS calculation |
| `from_proteins` | Predicted protein FASTA files | dbCAN HMM search, substrate mapping, and RPS calculation |
| `from_matrix` | Genome-by-feature CSV or TSV table with protein-family or substrate columns | Positive-value binarization, optional protein-family-to-substrate mapping, and RPS calculation |
| `download` | No biological input | Downloads the dbCAN HMM database used by FASTA workflows |

> **Important for FASTA modes:** Provide complete contigs/genomes or predicted protein FASTA files. Do not provide raw individual gene sequences unless they are intentionally being treated as a complete input dataset.

> **Important for matrix mode:** The first column is interpreted as the genome identifier by default. Every other column must be either a protein family or a substrate, as declared with `--matrix-type`. A table containing a list of substrates in one column is not the expected input format; features must be represented as columns.

---

## Matrix input mode

Matrix mode accepts a CSV or TSV table in which rows represent genomes and columns represent either **protein families** or **substrates**. The table must contain one genome identifier column and at least one feature column. A list-valued `substrates` column is not the input format for this mode.

Matrix mode bypasses gene prediction and HMM search. Substrate matrices do not require the substrate mapping, while protein-family matrices use data/substrate_key.json or the file supplied with -subs.

The input type must be declared explicitly with `--matrix-type`:

| Matrix type | Meaning of feature columns | Processing before RPS |
| :--- | :--- | :--- |
| `protein_families` | dbCAN/CAZy family identifiers such as `GH18`, `GH73`, `CE4`, or `AA1` | Positive family entries are binarized, then present families are mapped to substrates using `data/substrate_key.json`. |
| `substrates` | Substrate/resource names such as `chitin`, `peptidoglycan`, or `lignin` | Positive substrate entries are binarized and used directly as substrate presence features. |

Both binary values and counts are accepted. CaCo applies the following conversion before calculating RPS:

| Input value | Converted state |
| :--- | :--- |
| Greater than zero | Present, represented as `1` |
| Equal to zero | Absent, represented as `0` |
| Missing or nonnumeric | Absent, represented as `0` |

Counts are used only to determine presence. Their magnitude is not retained, so values of `1`, `10`, and `1,000` all represent the same present feature. In protein-family mode, multiple present families that map to the same substrate are collapsed into one substrate before the pairwise calculation.

### Protein-family matrix example

| genome | GH18 | GH73 | CE4 | AA1 |
| :--- | :--- | :--- | :--- | :--- |
| G1 | 2 | 0 | 1 | 0 |
| G2 | 1 | 3 | 0 | 1 |
| G3 | 0 | 0 | 0 | 0 |

In this example, `GH18` maps to host glycan, peptidoglycan, and chitin; `GH73` maps to peptidoglycan; `CE4` maps to xylan, peptidoglycan, and chitin; and `AA1` maps to lignin according to `data/substrate_key.json`. The resulting substrate sets are therefore derived from the mapping rather than from the protein-family names themselves.

Run it with:

```bash
python3 CaCo.py \
  -m from_matrix \
  -matrix tests/example/protein_family_matrix.tsv \
  --matrix-type protein_families \
  --matrix-id-column genome \
  -o carboncomp_output.tsv.xz \
  -cpus 8
```

### Substrate matrix example

| genome | chitin | peptidoglycan | lignin |
| :--- | :--- | :--- | :--- |
| G1 | 2 | 0 | 1 |
| G2 | 1 | 3 | 0 |
| G3 | 0 | 0 | 0 |

Run it with:

```bash
python3 CaCo.py \
  -m from_matrix \
  -matrix tests/example/substrate_matrix.tsv \
  --matrix-type substrates \
  --matrix-id-column genome \
  -o carboncomp_output.tsv.xz \
  -cpus 8
```

The `--matrix-id-column` option is optional. If omitted, CaCo uses the first column as the genome identifier. CSV and TSV delimiters are detected automatically. Matrix mode bypasses gene prediction and HMM search, but protein-family mode loads the substrate mapping from `data/substrate_key.json` unless a custom file is supplied with `-subs`.

---

## FASTA mode examples

### All-versus-all from nucleotide genomes

```bash
python3 CaCo.py \
  -m from_nucleotides \
  -gl example/gltest.txt \
  -cpus 16
```

The file supplied to `-gl` must contain one genome or protein path per line, depending on the selected mode.

### Pairwise comparison from nucleotide contigs

```bash
python3 CaCo.py \
  -m from_nucleotides \
  -g1 example/TARA_SAMEA2620259_METAG_HFKLHEHB.fa \
  -g2 example/TARA_SAMN05326651_METAG_RED00102.fa \
  -cpus 8
```

### Pairwise comparison from predicted proteins

```bash
python3 CaCo.py \
  -m from_proteins \
  -g1 example/TARA_SAMEA2620259_METAG_HFKLHEHB.faa \
  -g2 example/TARA_SAMN05326651_METAG_RED00102.faa \
  -cpus 8
```

### FASTA mode with Pyrodigal

```bash
python3 CaCo.py \
  -m from_nucleotides \
  -gl genomes.txt \
  -cpus 32 \
  --use-pyrodigal
```

---

## Command-line options

| Argument | Description |
| :--- | :--- |
| `-m` | Analysis mode: `download`, `from_proteins`, `from_nucleotides`, or `from_matrix`. |
| `-g1`, `-g2` | Paths to two genome or protein files for pairwise FASTA analysis. |
| `-gl` | File containing a list of genome or protein paths for all-versus-all FASTA analysis. |
| `-matrix` | Path to a genome-by-feature CSV or TSV table. Required when `-m from_matrix` is selected. |
| `--matrix-type` | Required for matrix mode. Select `protein_families` to map feature columns through the substrate key, or `substrates` to use substrate columns directly. |
| `--matrix-id-column` | Name of the column containing genome identifiers. Defaults to the first column. |
| `-o` | Output file name. Default: `carboncomp_output.tsv.xz`. The output is written to the current working directory unless an absolute path is supplied. |
| `-tmp` | Temporary directory for FASTA-mode intermediate files. Default: `tmp/`. |
| `-cpus` | Number of CPU cores to use. Default: all available cores. The workflow dynamically balances HMM-search threads to avoid oversubscription. |
| `--use-pyrodigal` | Use Pyrodigal for nucleotide FASTA gene prediction. If omitted, Prodigal is used. This option has no effect in matrix mode. |
| `-db` | Path to a custom dbCAN HMM database. Default: `data/dbcan.hmm`. Used by FASTA modes. |
| `-subs` | Path to a custom substrate-key JSON file. Default: `data/substrate_key.json`. Used by FASTA modes. |

---

## Output files

All final output files are written to the current working directory, meaning the directory from which CaCo is executed.

| File | Description |
| :--- | :--- |
| `allfams.tsv` | Feature table containing the present input features for each genome. In FASTA mode, these are dbCAN families. In matrix mode, these are the present protein-family or substrate columns. |
| `allsubs.tsv` | Substrate-list table used by the RPS calculation. In FASTA mode, these are inferred substrates. In matrix `protein_families` mode, these are produced by mapping present families through the substrate key. In matrix `substrates` mode, these are the present substrate columns. |
| `carboncomp_output.tsv.xz` | Compressed pairwise table containing overlap, competition, probability, RPS, and relative-RPS values. |

### Matrix-mode intermediate outputs

For the example matrix above, `allsubs.tsv` will have the following structure:

```text
genome\tsubstrates
G1\tacetate, xylose
G2\tacetate, glucose
G3\t
```

An empty feature list is valid. A genome with no positive annotations has an empty set of resources and is still included in pairwise comparisons.

### FASTA-mode feature outputs

`allfams.tsv` contains detected dbCAN families per genome:

| genome | families |
| :--- | :--- |
| TARA_SAMN05326651_METAG_RED00102 | AA3, PL12, GH102, CE11, GT2, GT4, GH73, GT19, AA7, GH109, GT5, GT41, CE4, GT9 |
| TARA_SAMEA2620259_METAG_HFKLHEHB | GH18, GT28, CE1, GT30, CE11, GT17, GH1, GT2, GT4, GT9, GH97, GH13, GT5, CE4, GH171, GH2 |

`allsubs.tsv` contains inferred substrates per genome:

| genome | substrates |
| :--- | :--- |
| TARA_SAMN05326651_METAG_RED00102 | cellooligosaccharide, host glycan, exo-polysaccharide, chitin, chitooligosaccharide, peptidoglycan, cellulose, xylan, glucooligosaccharide, lignin |
| TARA_SAMEA2620259_METAG_HFKLHEHB | alkaloid, alpha-glucan, alpha-mannan, arabinan, beta-fucosides, beta-galactan, beta-glucan, beta-glucuronan, beta-mannan, chitin, chitosan, exo-polysaccharide, glucosylglycerate, glycogen, host glycan, human milk polysaccharide, pectin, peptidoglycan, polyphenol, starch, sucrose, trehalose, uric acid, xylan |

### Pairwise competition output

The compressed pairwise output contains the following columns:

| Column | Description |
| :--- | :--- |
| `genome1`, `genome2` | Pair of genomes being compared. |
| `set1`, `set2` | Number of unique features detected in each genome. |
| `intersection` | Number of features shared by both genomes. |
| `competition` | Competition value calculated from the two feature sets and their intersection. |
| `relcomp` | Competition value normalized by the maximum possible competition for the pair. |
| `prob` | Estimated probability, based on 1,000 random overlap simulations, of observing at least the measured overlap under the implemented null model. |
| `RPS` | Resource Partitioning Score, calculated as `1 - 2 * competition`. |
| `relRPS` | Relative Resource Partitioning Score, calculated as `1 - 2 * relcomp`. |

Example row:

| genome1 | genome2 | set1 | set2 | intersection | competition | relcomp | prob | RPS | relRPS |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TARA_SAMN05326651_METAG_RED00102 | TARA_SAMEA2620259_METAG_HFKLHEHB | 10 | 24 | 5 | 0.1724 | 0.4138 | 0.349 | 0.6552 | 0.1724 |

---

## Interpretation of matrix-mode results

Matrix mode performs a **binary set-overlap analysis**. It does not use annotation abundance, gene copy number, expression level, or other count magnitude after conversion. In `protein_families` mode, RPS is calculated on the mapped substrate sets, not on the family names. If abundance-weighted comparisons are required, the current implementation would need a separate quantitative similarity or overlap model.

The probability calculation uses the implemented fixed feature-universe assumption in `CaCo.py`. Users applying matrix mode to a feature universe substantially different from the default assumption should validate the probability interpretation for their dataset.

---

## Performance

The FASTA workflow is fully parallelized. For thousands of genomes, wall time can decrease substantially on a modern multi-core server.

| Step | Original | Current implementation |
| :--- | :--- | :--- |
| Gene prediction | Sequential Prodigal | Parallel Prodigal or Pyrodigal |
| HMM search | Sequential | Parallel with dynamically balanced threads |
| Parsing | Sequential | Parallel |
| Feature extraction | Sequential | Parallel |
| Pairwise RPS | Sequential | Chunked and parallel |
| Matrix conversion | Not available | Local CSV/TSV loading followed by binary conversion |

See [What's New](whats_new.md) for additional performance details.

---

## Citation

If you use CaCo in your research, please cite:

> Santos-Júnior, C. et al. (2026). Resource availability structures microbial competition through genomic niche partitioning. *PNAS*, **123**(18): e2526391123.

---

## License

This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

---

## Issues and support

Please report bugs or feature requests through the [GitHub issue tracker](https://github.com/celiosantosjr/CaCo/issues).
