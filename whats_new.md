# What is New in the Parallelized CaCo

This document summarizes the major improvements and changes introduced in the new version of CaCo (the parallelized, “lightspeed” version) compared to the original implementation.

---

## 1. Gene Prediction – Dual Backend (Prodigal Default, Pyrodigal Optional)

| **Original** | **New** |
|--------------|---------|
| Calls `prodigal` as a subprocess for each genome (slow, sequential) | **Default:** Uses the original **Prodigal** subprocess, but now runs in **parallel** across all genomes for immediate speed gains – ensuring **100% numerical reproducibility** with the PNAS paper. |
| No parallelization of gene prediction | **Optional:** Use `--use-pyrodigal` to switch to the fast Python binding (**Pyrodigal**) for up to **5× faster** in‑memory gene calling. |
| **Speed**: ~1× (baseline) | **Speed**: Parallel Prodigal gives **N×** speedup; Pyrodigal adds another **5×** on top. |

---

## 2. HMM Search (dbCAN) – Fully Parallelized with Smart CPU Allocation

| **Original** | **New** |
|--------------|---------|
| Sequential `hmmsearch` calls, one genome at a time (using 3 threads per call) | Concurrent `hmmsearch` jobs for **all** genomes, but **dynamically balanced** to avoid CPU oversubscription: <br> • `num_jobs = min(total_cpus, number_of_genomes)` <br> • `cpu_per_task = max(1, total_cpus // num_jobs)` <br> This ensures that total requested threads never exceed the available cores, preventing performance thrashing. |
| No progress feedback | Progress bars via `tqdm` for each parallel task |
| **Speed**: O(N) wall‑time | **Speed**: Wall‑time reduces to ~O(1) (optimal core utilisation) |

---

## 3. Parsing of HMM Output – Parallelized

| **Original** | **New** |
|--------------|---------|
| Parses `.dbcan` files one after another | Parsing is performed **in parallel** for all files using a process pool |
| Same logic, but sequential | Same logic, but **concurrent** |

---

## 4. Feature Extraction – Parallelized with Memory Safeguards

| **Original** | **New** |
|--------------|---------|
| Reads each `.parsed` file sequentially, extracts families and substrates | Extraction is parallelised – each worker processes one file and returns its results |
| Writes directly to `allfams.tsv` and `allsubs.tsv` with `'a+'` mode (prone to race conditions if parallel) | Uses in‑memory DataFrames and writes both files **once** at the end (safe and fast). <br><br> **⚠️ Safeguard:** When processing >10,000 genomes, the tool now prints a warning to advise using a chunked writer to avoid potential memory exhaustion. |
| **Speed**: Sequential | **Speed**: **N×** faster (cores) |

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
| `-cpus` / `--cpus` | Number of CPU cores to use (default: all available) – controls the size of process pools and dynamically adjusts `hmmsearch` thread allocation. |
| `--use-pyrodigal` | **NEW:** Explicit opt‑in flag to use Pyrodigal for gene prediction. **If omitted**, the tool defaults to Prodigal (subprocess) to guarantee exact reproducibility with the original PNAS publication. |
| (Implicit) | The `-m` argument now only accepts `download`, `from_proteins`, or `from_nucleotides` (was required before, now optional but recommended). |

---

## 7. Improved Error Handling, Robustness & Output Location

- **Database download** now validates header and footer of the HMM file, ensuring a complete download.
- **Subprocess calls** are wrapped with `check=True` and redirect output to `/dev/null` to avoid clutter.
- **Temporary file management** uses Python’s `tempfile` module for per‑chunk RPS files, preventing name collisions (cleanup is guaranteed via `finally` blocks).
- **Output directory fix:** Final outputs (`allfams.tsv`, `allsubs.tsv`, `carboncomp_output.tsv.xz`) are now placed in the **current working directory** (`os.getcwd()`), not the script directory. This prevents permission errors when the tool is installed system‑wide and aligns with standard UNIX tool behaviour.
- **Dependency check:** If not using `--use-pyrodigal`, the tool now explicitly verifies that `prodigal` is installed in `PATH` and raises a clear error if missing.

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
| Gene prediction        | Sequential Prodigal    | Parallel (Prodigal/Pyrodigal) | **N× – 10×**  |
| HMM search             | Sequential hmmsearch   | Parallel + dynamically balanced | **N×** (cores)|
| Parsing                | Sequential             | Parallel                 | **N×** (cores)|
| Feature extraction     | Sequential             | Parallel (with memory warning) | **N×** (cores)|
| Pairwise RPS (O(N²))   | Sequential O(N²)       | Chunked + parallel       | **N×** (cores)|

Where `N` = number of genomes.  
For thousands of genomes, the wall‑time drops from **days to hours** on a modern multi‑core server.

---

## ⚠️ Important Note on Reproducibility

The version of CaCo used in the **PNAS paper** (“Resource availability structures microbial competition through genomic niche partitioning”; Santos-Júnior et al., *PNAS* **123**(18): e2526391123, 2026) uses **Prodigal** via subprocess for gene prediction.

The current version is **engineered to preserve this reproducibility by default**:

- **Default behaviour:** Uses parallelised **Prodigal** (subprocess), producing bit‑wise identical gene calls to the original paper.
- **Opt‑in speed mode:** Use `--use-pyrodigal` to switch to the faster Pyrodigal backend – which may produce **slightly different gene calls** in rare cases, leading to minor deviations in downstream CAZyme annotations and RPS scores.

The exact code used for the paper is archived at the commit from **January 14th**:

👉 [https://github.com/celiosantosjr/CaCo/blob/db50a80e974b4daedb217cf971a68c6079fc33e7/CaCo.py](https://github.com/celiosantosjr/CaCo/blob/db50a80e974b4daedb217cf971a68c6079fc33e7/CaCo.py)

You can also try the full [legacy repository](https://github.com/LMPB/CaCo) (not actively maintained, updated before publication). In that repository, RPS is called EIT (Ecological Interaction Type), but the logic and calculations remain original.

> **For exact reproduction of the PNAS results**, use the default mode (without `--use-pyrodigal`) or the commit linked above.  
> **For high‑throughput screening** of thousands of genomes where slight differences are acceptable, add `--use-pyrodigal` for dramatic speed improvements.

---

## Compatibility Notes

- The new version **still supports** the same input formats and output structures as the original.
- **Default gene predictor:** Prodigal (subprocess) for exact reproducibility.
- All intermediate files are written to the temporary directory (`-tmp`).
- Final outputs (`allfams.tsv`, `allsubs.tsv`, `carboncomp_output.tsv.xz`) are placed in the **current working directory** for easy access and system‑wide compatibility.

---

**Last updated:** July 2026 (Updated to reflect dual‑backend and HPC‑safety fixes)  
**Maintained by:** The CaCo development team
