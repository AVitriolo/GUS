# METEOR on TCGA, one cancer type at a time (gene level, SLURM)

The analysis done for TCGA PRAD (`/home/ieo7293/TCGA_PRAD`), made generic per TCGA cohort. The model settings are those of the prostate WGBS runs.

| Step | What | Same as TCGA PRAD? |
|---|---|---|
| 0 | Select cohorts on the UCSC Xena tcgaHub with **≥ 100 samples having both** 450K methylation and HiSeqV2 expression. Combined cohorts (PANCAN, COADREAD, GBMLGG, LUNG, …) are excluded. | new |
| 1 | Download `HumanMethylation450` and `HiSeqV2` (gene-level RSEM, log2(x+1)) from Xena | yes (`code/0_download_data.R`) |
| 2 | `data.processing()`: beta → M (0/1 → −10/10), REMP `grooMethy` imputation of probes with ≤ 20% NA, drop probes still NA, keep genes with mean expression > 1, matched samples, hg19 annotation, NCBI synonym search; export CSVs and gene list | yes (unchanged logic, only paths are arguments) |
| 2b | After preprocessing, skip the cohort if fewer than 100 samples remain matched | new |
| 3 | Per gene: TSS anchor = probe with the smallest \|distToTSS\| (Zhou-lab hg19 gene association); features = all probes within **±1 Mb**; target = gene expression; all matched samples, **70/30 split** with seed 48047510 | yes (`make_datasets_prad_positional_cli.py`; tables verified identical on PRAD) |
| 4 | SVR + SelectKBest, **RepeatedKFold 3 × 5, k capped at the per-fold n**, RandomSearchCV 27 iterations, no SHAP | **as the prostate WGBS runs** (`model_selection_split_kbest_repeated_kcapped.py`); the original PRAD run used single 3-fold and uncapped k |

## Setup (once)

1. Copy this folder to the cluster, e.g. `/beegfs/scratch/ieo7293/tcga_meteor_pancancer`.
2. Create the environment:
   ```
   conda env create -f env/meteor_tcga.yml
   ```
   Or set `CONDA_ENV` in the config to the path of an existing environment.
3. Edit `config/pipeline.conf`:
   - `WORK_DIR`: where data and results go, about 10–20 GB per cohort;
   - `SLURM_PARTITION` / `SLURM_ACCOUNT` / `SLURM_QOS`: leave empty to use the cluster's defaults;
   - the resource settings, if needed.

## Run

Everything is submitted with `sbatch`; no work happens on the login node. From the pipeline folder:

```
mkdir -p logs
sbatch submit_all.sbatch                  # selects the cohorts and submits all of them (1 CPU, 4 GB, about 10-30 min)
tail -f logs/submit_all_<jobid>.out       # progress: one line per step and per cohort while it counts the samples
./status.sh                               # any time: per cohort the step (queue / download / preprocessing / genes / finished), genes done / total, failures
```

- To submit chosen cohorts only, use `./submit_cohort.sh BRCA LUAD`. It only calls `sbatch`, so it's light.
- The submit job first checks that the compute node can reach the Xena server, and stops with a clear message if it can't. The downloads need internet.
- Each cohort runs as one **prepare** job (download + preprocessing, 8 CPUs, 128 GB). At its end, the prepare job submits a **gene job array**: 100 genes per task, 8 models at a time × 2 cores, 48 GB, at most 20 tasks at once per cohort. All cohorts run in parallel, as far as your cluster quota allows.
- Rerunning is safe: finished downloads, preprocessing and genes are skipped.
- Logs: `$WORK_DIR/<COHORT>/logs/prepare_*.out` and `genes_*_*.out`. `python scripts/collect_results.py $WORK_DIR BRCA LUAD` writes the per-gene summary tables.

## Outputs (per cohort, in `$WORK_DIR/<COHORT>/`)

| Path | Content |
|---|---|
| `results/performance/<GENE>_SVR.tsv` | train / cv / test R² and other metrics |
| `results/sel_feats/<GENE>_SVR.tsv` | selected probes |
| `results/hyper/<GENE>_SVR.tsv` | best hyperparameters |
| `summary_<COHORT>.tsv` | one row per gene (written by `collect_results.py`) |
| `<COHORT>_data_meth.csv`, `need_gene_exp_<COHORT>.csv`, … | preprocessed data, same format as the TCGA PRAD files |
| `logs/` | prepare and array logs; `results/logs/<GENE>.log` per model |

Model `.pkl` files and the per-gene train/test tables are deleted after each gene, to save space. Set `KEEP_GENE_DATASETS=1` to keep the tables.

## Notes

- Run time: about 5–10 minutes per gene on 2 cores, measured on PRAD (405 fits per model with 3 × 5 repeated CV and 27 iterations). That's about 1.5–2.5 h per array task of 100 genes. A cohort with about 14,000 genes is 140 tasks, about 10–18 h with 20 tasks in parallel, plus 2–6 h of preprocessing: about 1 day per cohort. All cohorts together need about 80–100k core-hours, about 11 days at 320 cores or about 4 days at 1,000. Raise `MAX_CONCURRENT_TASKS` if your quota allows; `N_REPEATS=3` roughly halves the time.
- The REMP imputation in step 2 is the slowest and most memory-hungry part for large cohorts (BRCA, about 800 samples). Raise `PREP_MEM` / `PREP_TIME` if the prepare job is killed.
- Like the PRAD analysis, all matched samples are used (tumour and normal), coordinates are hg19 (Xena), and gene names come from Xena after the synonym search.

## Stability selection (after the gene models)

The same procedure as the prostate WGBS analysis, per cohort. For each gene model:
- **Rebuild the train table** with `03_make_gene_datasets.py`: same construction and seed, so the same train samples and candidate probes.
- **Half-samples:** `STAB_N_PAIRS` (100) complementary pairs of the train samples, each half stratified by sample type from the barcode (tumour / normal).
- **Selection on each half:** `SelectKBest(f_regression)` with the model's k.
- **Per probe:** selection frequency, and stable if it passes the `stabs` cutoff for PFER ≤ `STAB_PFER`.
- **Per gene:** Nogueira stability Φ (`stabm`).
- **Per co-methylated region** (adjacent r ≥ `STAB_CORR`, ≤ `STAB_MAX_GAP` bp, hg19 probe positions): the same, at region level.

Setup, once: `conda env create -f env/meteor_stability.yml`, or set `STAB_R_ENV` in the config to an existing env with R, stabm, stabs and data.table.

```
./submit_stability.sh BRCA LUAD            # chosen cohorts (gene models must be finished; --force to override)
./submit_stability.sh all                  # every cohort whose gene models are finished
./stability_status.sh                      # chunks done / total, jobs, summary written
./submit_stability.sh --summary_only PRAD  # only the summary, e.g. after changing STAB_PFER
```

- Each cohort runs one array task per gene chunk (2 CPUs, 48 GB): it rebuilds the chunk's tables and runs the selections. A summary job then runs the two R steps and the report.
- Finished chunks are skipped on rerun.
- **PFER and array data.** With about 400 candidate probes and k about 100, the cutoff for PFER ≤ 1 is about 1.0, so few probes can be called stable however reproducible they are. Frequencies and Φ are always written, so `STAB_PFER` can be changed afterwards by rerunning only the summary.

Outputs in `$WORK_DIR/<COHORT>/stability_<scale>/`:

| File | Content |
|---|---|
| `stability_report.png` / `.pdf` | % of each model's selected probes in a stable region, and Φ, by CV R² bin |
| `stability_report_genes.tsv` | one row per gene: cv / test R², k, candidates, Φ, stable probes and regions |
| `cpg_freq.tsv` | per probe: selection frequency, saved_selected, stable |
| `region_freq_corr<r>.tsv`, `region_genes_corr<r>.tsv` | per region and per gene, region level |
| `stability_genes.tsv` | per gene: Φ, PFER cutoff, number of stable probes |
