#!/usr/bin/env bash
# Manual run: bMIND (reference-free) purification of CPGEA Tumour RNA with
# PAMES purity. One joint fit on all Tumours (train + test together, as for the
# metastases); all Tumours written. Not wired into any Snakefile - run by hand.
set -euo pipefail
R_ENV=/home/ieo7293/miniforge3/envs/r_mind/bin/Rscript
NCORE=${NCORE:-16}
D=results/bmind/CPGEA
mkdir -p $D logs/bmind
$R_ENV analyses/scripts/run_bmind.R \
    --counts=resources/RNA/kallisto_counts_NT.tsv \
    --transcripts=resources/RNA/kallisto_counts_NTM_train_25_10.tsv \
    --purity=resources/purity/pames_purity_NT.tsv \
    --fit_samples=$D/all_tumours.txt \
    --out_prefix=$D/all \
    --ncore=$NCORE > logs/bmind/CPGEA_all.log 2>&1
echo "[$(date)] done"
