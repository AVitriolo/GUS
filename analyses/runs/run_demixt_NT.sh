#!/usr/bin/env bash
# Manual run: DeMixT purification of CPGEA Tumour RNA with the setting chosen
# on the train grid (run_demixt_grid_NT.sh): 10 spike-ins, 1500 genes.
# Train proportions are reused from the grid cache; test is estimated fresh
# (train Normals as reference for both). Counts +1 before DeMixT_S2.
# Not wired into any Snakefile - run by hand.
set -euo pipefail
R_ENV=/home/ieo7293/miniforge3/envs/r_demixt/bin/Rscript
THREADS=${THREADS:-16}
mkdir -p results/demixt_NT/final logs/demixt
[ -f results/demixt_NT/final/demixt_GS_train_ns10_ng1500.rds ] || \
    cp results/demixt_NT/grid/demixt_GS_train_ns10_ng1500.rds results/demixt_NT/final/
OMP_NUM_THREADS=$THREADS $R_ENV analyses/scripts/run_demixt_NT.R \
    --counts=resources/RNA/kallisto_counts_NT.tsv \
    --transcripts=resources/RNA/kallisto_counts_NTM_train_25_10.tsv \
    --coldata=resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds \
    --outdir=results/demixt_NT/final --splits=train,test --stage=all \
    --nspikein=10 --ngene=1500 --nthread=$THREADS
