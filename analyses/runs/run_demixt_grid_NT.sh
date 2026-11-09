#!/usr/bin/env bash
# Manual run: DeMixT parameter grid on the TRAIN split (proportion step only),
# as recommended in the DeMixT vignette: spike-ins x number of selected genes.
# The setting is chosen on train only, then applied to train and test
# (run_demixt_NT.R --stage all).
#
# Spike-ins capped at min(n/3, 200) = 29 for 88 train Tumours.
# Not wired into any Snakefile - run by hand. Cached runs are reused.
#
# Usage: ./run_demixt_grid_NT.sh            (N_JOBS=4 x THREADS=6)
set -euo pipefail

R_ENV=/home/ieo7293/miniforge3/envs/r_demixt/bin/Rscript
N_JOBS=${N_JOBS:-4}
THREADS=${THREADS:-6}
OUTDIR=results/demixt_NT/grid
LOGDIR=logs/demixt/grid
mkdir -p $OUTDIR $LOGDIR

run_one() {
    ns=$1; ng=$2
    OMP_NUM_THREADS=$THREADS $R_ENV analyses/scripts/run_demixt_NT.R \
        --counts=resources/RNA/kallisto_counts_NT.tsv \
        --transcripts=resources/RNA/kallisto_counts_NTM_train_25_10.tsv \
        --coldata=resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds \
        --outdir=$OUTDIR --splits=train --stage=gs \
        --nspikein=$ns --ngene=$ng --nthread=$THREADS \
        > $LOGDIR/train_ns${ns}_ng${ng}.log 2>&1 \
        && echo "[$(date +%H:%M)] ns=$ns ng=$ng done" \
        || echo "[$(date +%H:%M)] ns=$ns ng=$ng FAILED (see $LOGDIR/train_ns${ns}_ng${ng}.log)"
}
export -f run_one
export R_ENV THREADS OUTDIR LOGDIR

for ns in 0 10 20 29; do for ng in 500 1000 1500 2500; do echo "$ns $ng"; done; done \
    | xargs -P $N_JOBS -L 1 bash -c 'run_one $0 $1'

echo "[$(date)] grid finished: $(ls $OUTDIR/demixt_pi_train_*.tsv | wc -l)/16 settings"
