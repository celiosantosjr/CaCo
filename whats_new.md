# What is New in the Parallelized CaCo

This document summarizes the major improvements and changes introduced in the new version of CaCo (the parallelized, “lightspeed” version) compared to the original implementation.

---

## 1. Gene Prediction – Replaced Prodigal with Pyrodigal

| **Original** | **New** |
|--------------|---------|
| Calls `prodigal` as a subprocess for each genome (slow, sequential) | Uses **Pyrodigal** – a fast Python binding of Prodigal – called directly in‑memory |
| No parallelization of gene prediction | Gene prediction runs in **parallel** across all genomes using `ProcessPoolExecutor` |
| **Speed**: ~1× (baseline) | **Speed**: up to **5× faster** for gene calling (CPU‑bound) |

---

## 2. HMM Search (dbCAN) – Fully Parallelized

| **Original** | **New** |
|--------------|---------|
| Sequential `hmmsearch` calls, one genome at a time (using only 3 threads per call) | Concurrent `hmmsearch` jobs for **all** genomes, each using `--cpu 2` (configurable) |
| No progress feedback | Progress bars via `tqdm` for each parallel task |
| **Speed**: O(N) wall‑time | **Speed**: Wall‑time reduces to ~O(1) (assuming enough cores) |

---

## 3. Parsing of HMM Output – Parallelized

| **Original** | **New** |
|--------------|---------|
| Parses `.dbcan` files one after another | Parsing is performed **in parallel** for all files using a process pool |
| Same logic, but sequential | Same logic, but **concurrent** |

---

## 4. Feature Extraction – Parallelized

| **Original** | **New** |
|--------------|---------|
| Reads each `.parsed` file sequentially, extracts families and substrates | Extraction is parallelised – each worker processes one file and returns its results |
| Writes directly to `allfams.tsv` and `allsubs.tsv` with `'a+'` mode (prone to race conditions if parallel) | Uses in‑memory DataFrames, then writes both files **once** at the end (safe and fast) |

---

## 5. Pairwise RPS Computation – Chunked and Parallel

| **Original** | **New** |
|--------------|---------|
| Brute‑force O(N²) **sequential** loop over all genome pairs (slow for >1000 genomes) | Pairs are split into **chunks** and processed in **parallel** by multiple workers |
| Writes directly to the final compressed file (locks / IO contention) | Each worker writes to a **temporary file**; final merge is done after all chunks complete – minimal IO contention |
| **Speed**: O(N²) wall‑time | **Speed**: Wall‑time reduced by factor = number of workers (ideal scaling) |

---

## 6. Added Command‑Line Options

| **New Argument** | **Description** |
|------------------|-----------------|
| `-cpus` / `--cpus` | Number of CPU cores to use (default: all available) – controls the size of process pools |
| (Implicit) | The `-m` argument now only accepts `download`, `from_proteins`, or `from_nucleotides` (was required before, now optional but recommended) |

---

## 7. Improved Error Handling & Robustness

- **Database download** now validates header and footer of the HMM file, ensuring a complete download.
- **Subprocess calls** are wrapped with `check=True` and redirect output to `/dev/null` to avoid clutter.
- **Temporary file management** uses Python’s `tempfile` module for per‑chunk RPS files, preventing name collisions.
- **File paths** are resolved relative to the script directory (`SCRIPT_DIR`), ensuring `allfams.tsv` and `allsubs.tsv` are written in the expected location for testing.

---

## 8. Code Organization & Maintainability

- Clear separation of concerns: each parallel task has its own function (e.g., `parallel_hmmsearch`, `parallel_parse`, `extract_features_parallel`).
- Consistent use of `ProcessPoolExecutor` across all parallel steps.
- Use of `tqdm` to provide real‑time progress feedback for each major stage.
- Modular structure makes it easy to swap in alternative predictors or HMM search engines.

---

## 9. GitHub Actions Integration

- The repository now includes **CI/CD workflows** (`.github/workflows/`).
- Automated tests are run on **every push** and **pull request**.
- Tests validate:
  - Gene prediction modes (`from_nucleotides`, `from_proteins`)
  - Pairwise and list‑based input
  - Output files (`allfams.tsv`, `allsubs.tsv`, `carboncomp_output.tsv.xz`) against expected results.
- This ensures that any future changes do not break the core functionality.

---

## Summary of Performance Gains

| **Step**               | **Original**           | **Parallelized**         | **Speedup**    |
|------------------------|------------------------|--------------------------|----------------|
| Gene prediction        | Sequential Prodigal    | Parallel Pyrodigal       | **5× – 10×**  |
| HMM search             | Sequential hmmsearch   | Parallel hmmsearch       | **N×** (cores)|
| Parsing                | Sequential             | Parallel                 | **N×** (cores)|
| Feature extraction     | Sequential             | Parallel                 | **N×** (cores)|
| Pairwise RPS (O(N²))   | Sequential O(N²)       | Chunked + parallel       | **N×** (cores)|

Where `N` = number of genomes.  
For thousands of genomes, the wall‑time drops from **days to hours** on a modern multi‑core server.

---

## ⚠️ Important Note on Reproducibility

The version of CaCo used in the **PNAS paper** (“Resource availability structures microbial competition through genomic niche partitioning”; Santos-Júnior et al., *PNAS* **123**(18): e2526391123, 2026) is **not 100% identical** to the current parallelized version.

The exact code used for the paper is archived at the commit from **January 14th**:
👉 [https://github.com/celiosantosjr/CaCo/blob/db50a80e974b4daedb217cf971a68c6079fc33e7/CaCo.py](https://github.com/celiosantosjr/CaCo/blob/db50a80e974b4daedb217cf971a68c6079fc33e7/CaCo.py)

**Why the difference?**  
The original paper version used **Prodigal** for gene prediction via subprocess calls. The current version replaces Prodigal with **Pyrodigal** (a fast Python binding) to achieve parallelization and speed gains. While Pyrodigal is functionally similar, it may produce **slightly different gene calls** than Prodigal in some cases, which can lead to minor deviations in downstream CAZyme annotations and RPS scores.

> **For exact reproduction of the PNAS results**, please use the commit linked above.  
> **For high‑throughput screening** of thousands of genomes where slight differences are acceptable, the current parallelized version offers dramatic speed improvements.

---

## Compatibility Notes

- The new version **still supports** the same input formats and output structures as the original.
- The default gene predictor is now **Pyrodigal**.
- All intermediate files are now written to the temporary directory (`-tmp`), while the final outputs (`allfams.tsv`, `allsubs.tsv`, `carboncomp_output.tsv.xz`) are placed in the script’s directory for easy access and testing.

---

**Last updated:** July 2026  
**Maintained by:** The CaCo development team
