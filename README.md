# CaCo

[![DBcan Database](https://img.shields.io/badge/DBcan-Database-blue)](https://bcb.unl.edu/dbCAN2/)
[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub issues](https://img.shields.io/badge/Issues-Open-red)](https://github.com/celiosantosjr/CaCo/issues)
[![CI](https://github.com/celiosantosjr/CaCo/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/celiosantosjr/CaCo/actions)

**CaCo**: A high‑throughput program for predicting carbon source competition and ecological interaction types from microbial genomes.

> **New in this version:** Fully parallelized workflow, smart CPU allocation, and **dual‑backend gene prediction** – defaulting to Prodigal for exact reproducibility with the PNAS paper, with an opt‑in Pyrodigal mode for maximum speed.  
> 👉 See the [What's New](whats_new.md) document for all performance gains and changes.

---

## 📌 Important Note on Reproducibility

The version of CaCo used in the **PNAS paper** (“Resource availability structures microbial competition through genomic niche partitioning”; Santos-Júnior et al., *PNAS* **123**(18): e2526391123, 2026) uses **Prodigal** via subprocess for gene prediction.

- **Default behaviour (now):** Uses **Prodigal** (subprocess) **in parallel**, guaranteeing bit‑wise identical gene calls to the original publication.
- **Opt‑in speed mode:** Add `--use-pyrodigal` to use Pyrodigal (fast Python binding). This is **up to 5× faster** but may produce slightly different gene calls in rare cases, leading to minor deviations in downstream scores.

> For exact reproduction of the PNAS results, run without `--use-pyrodigal` or use the [archived commit](https://github.com/celiosantosjr/CaCo/blob/db50a80e974b4daedb217cf971a68c6079fc33e7/CaCo.py) from January 14th.

---

## Installation

First, clone the repository:

```bash
git clone https://github.com/celiosantosjr/CaCo
cd CaCo/
```

Set up the Conda environment (this installs all dependencies, including **Prodigal** and **Pyrodigal**):

```bash
conda env create -f CaCo_env.yml
conda activate CaCo
```

> If you prefer to install Prodigal manually: `conda install -c bioconda prodigal`

Download the dbCAN HMM database (required for CAZyme annotation):

```bash
python3 CaCo.py -m download
```

Check that everything works:

```bash
python3 CaCo.py -h
```

---

## Usage

CaCo accepts either:
- Two genomes (pairwise mode), or
- A file listing multiple genome paths (all‑vs‑all mode).

It works with both **nucleotide contigs** (`from_nucleotides`) and **predicted protein files** (`from_proteins`).

> **⚠️ Important:** Provide complete contigs/genomes or protein FASTA files – **do not** input raw gene sequences.

### Basic examples

**1. All‑vs‑all from a list of nucleotide genomes:**

```bash
python3 CaCo.py -m from_nucleotides -gl example/gltest.txt -cpus 16
```

**2. Pairwise comparison from nucleotide contigs:**

```bash
python3 CaCo.py -m from_nucleotides \
  -g1 example/TARA_SAMEA2620259_METAG_HFKLHEHB.fa \
  -g2 example/TARA_SAMN05326651_METAG_RED00102.fa \
  -cpus 8
```

**3. Pairwise comparison from predicted proteins (skip gene calling):**

```bash
python3 CaCo.py -m from_proteins \
  -g1 example/TARA_SAMEA2620259_METAG_HFKLHEHB.faa \
  -g2 example/TARA_SAMN05326651_METAG_RED00102.faa \
  -cpus 8
```

**4. Opt‑in for maximum speed (Pyrodigal):**

```bash
python3 CaCo.py -m from_nucleotides -gl genomes.txt -cpus 32 --use-pyrodigal
```

---

## Command‑Line Options

| Argument | Description |
| :--- | :--- |
| `-m` | Mode: `download`, `from_proteins`, or `from_nucleotides`. |
| `-g1`, `-g2` | Paths to two genomes (pairwise mode). |
| `-gl` | File containing a list of genome paths (all‑vs‑all mode). |
| `-o` | Output file name (default: `carboncomp_output.tsv.xz`). **Placed in the current working directory.** |
| `-tmp` | Temporary directory for intermediate files (default: `tmp/`). |
| `-cpus` / `--cpus` | Number of CPU cores to use (default: all available). The tool dynamically balances threads to avoid oversubscription. |
| `--use-pyrodigal` | **Optional:** Use Pyrodigal for gene prediction (faster, but slightly different from original Prodigal). If omitted, the script uses Prodigal (subprocess) for exact reproducibility. |
| `-db` | Path to custom dbCAN HMM database (default: `data/dbcan.hmm`). |
| `-subs` | Path to substrate key JSON (default: `data/substrate_key.json`). |

---

## Output Files

All final output files are written to the **current working directory** (the directory from which you run the script):

| File | Description |
| :--- | :--- |
| `allfams.tsv` | Table of dbCAN protein families detected in each genome. |
| `allsubs.tsv` | Table of carbon substrates predicted to be consumed by each genome. |
| `carboncomp_output.tsv.xz` | Compressed table of pairwise competition scores (RPS, relative RPS, probabilities, etc.). |

### Output examples

**`allfams.tsv`** – families per genome:

| genome | families |
| :--- | :--- |
| TARA_SAMN05326651_METAG_RED00102 | AA3, PL12, GH102, CE11, GT2, GT4, GH73, GT19, AA7, GH109, GT5, GT41, CE4, GT9 |
| TARA_SAMEA2620259_METAG_HFKLHEHB | GH18, GT28, CE1, GT30, CE11, GT17, GH1, GT2, GT4, GT9, GH97, GH13, GT5, CE4, GH171, GH2 |

**`allsubs.tsv`** – substrates per genome:

| genome | substrates |
| :--- | :--- |
| TARA_SAMN05326651_METAG_RED00102 | cellooligosaccharide, host glycan, exo-polysaccharide, chitin, chitooligosaccharide, peptidoglycan, cellulose, xylan, glucooligosaccharide, lignin |
| TARA_SAMEA2620259_METAG_HFKLHEHB | alkaloid, alpha-glucan, alpha-mannan, arabinan, beta-fucosides, beta-galactan, beta-glucan, beta-glucuronan, beta-mannan, chitin, chitosan, exo-polysaccharide, glucosylglycerate, glycogen, host glycan, human milk polysaccharide, pectin, peptidoglycan, polyphenol, starch, sucrose, trehalose, uric acid, xylan |

**`carboncomp_output.tsv.xz`** – pairwise competition scores (columns):

| Column | Description |
| :--- | :--- |
| `genome1`, `genome2` | Pair of genomes. |
| `set1`, `set2` | Number of unique carbon substrates detected in each. |
| `intersection` | Number of overlapping substrates. |
| `competition` | Jaccard distance between the two substrate sets. |
| `relcomp` | Jaccard distance normalized by the maximum possible distance. |
| `prob` | Permutation probability (1,000 simulations) of observing such overlap by chance. |
| `RPS` | Resource Partitioning Score = `1 - 2 * competition`. |
| `relRPS` | Relative RPS = `1 - 2 * relcomp`. |

Example row:

| genome1 | genome2 | set1 | set2 | intersection | competition | relcomp | prob | RPS | relRPS |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| TARA_SAMN05326651_METAG_RED00102 | TARA_SAMEA2620259_METAG_HFKLHEHB | 10 | 24 | 5 | 0.1724 | 0.4138 | 0.349 | 0.6552 | 0.1724 |

---

## Performance

The workflow is fully parallelized. For thousands of genomes, wall‑time drops from **days to hours** on a modern multi‑core server.

| Step | Original | Parallelized | Speedup |
| :--- | :--- | :--- | :--- |
| Gene prediction | Sequential Prodigal | Parallel (Prodigal or Pyrodigal) | **N× – 10×** |
| HMM search | Sequential | Parallel + dynamically balanced threads | **N×** (cores) |
| Parsing | Sequential | Parallel | **N×** (cores) |
| Feature extraction | Sequential | Parallel | **N×** (cores) |
| Pairwise RPS (O(N²)) | Sequential | Chunked + parallel | **N×** (cores) |

See the [What's New](whats_new.md) document for full details.

---

## Citation

If you use CaCo in your research, please cite:

> Santos-Júnior, C. et al. (2026). Resource availability structures microbial competition through genomic niche partitioning. *PNAS*, **123**(18): e2526391123.

---

## License

This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

---

## Issues & Support

Please report bugs or feature requests via the [GitHub issue tracker](https://github.com/celiosantosjr/CaCo/issues).
