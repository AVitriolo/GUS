#!/usr/bin/env bash
# Manual run: bMIND (reference-free) purification of WCDT metastatic RNA with
# PAMES purity (pames_purity_M.tsv), one fit per biopsy site (liver excluded).
# Not wired into any Snakefile - run by hand.
set -euo pipefail
R_ENV=/home/ieo7293/miniforge3/envs/r_mind/bin/Rscript
NCORE=${NCORE:-16}
D=results/bmind/M
mkdir -p $D logs/bmind
for SITE in Bone Lymph_node Other; do
    $R_ENV analyses/scripts/run_bmind.R \
        --counts=resources/RNA/kallisto_counts_M.tsv \
        --transcripts=resources/RNA/kallisto_counts_NTM_train_25_10.tsv \
        --purity=resources/purity/pames_purity_M.tsv \
        --fit_samples=$D/${SITE}_samples.txt \
        --out_prefix=$D/$SITE \
        --ncore=$NCORE > logs/bmind/M_${SITE}.log 2>&1
    echo "[$(date +%H:%M)] $SITE done: $(grep '^written' logs/bmind/M_${SITE}.log)"
done
