#!/usr/bin/env bash
# Stability selection, one array task = one gene chunk (chunks/chunk_<task id>.txt, the chunks of the gene array).
# Submitted by submit_stability.sh. A finished chunk (stability_<scale>/done_chunk_<id>) is skipped when rerun.
set -uo pipefail
COHORT=$1
source "$PIPE_DIR/config/pipeline.conf"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
sleep $((RANDOM % 60))   # spread the start of the array tasks (BeeGFS)
C=$WORK_DIR/$COHORT
CHUNK=$(printf "%s/chunks/chunk_%04d.txt" "$C" "$SLURM_ARRAY_TASK_ID")
[ -f "$CHUNK" ] || { echo "no chunk $CHUNK"; exit 1; }
LOG=$C/logs/stability_chunk_$(printf "%04d" "$SLURM_ARRAY_TASK_ID").log
# retry when python fails while importing packages (BeeGFS hiccup), not on real errors (as run_gene_chunk.sh)
for try in 1 2 3; do
  python "$PIPE_DIR/scripts/stability_gene_chunk.py" --cohort=$COHORT --work_dir=$C --pipe_dir=$PIPE_DIR --chunk=$CHUNK \
    --scale=$METH_SCALE --window_bp=$WINDOW_BP --test_size=$TEST_SIZE --n_pairs=$STAB_N_PAIRS --corr=$STAB_CORR \
    --max_gap=$STAB_MAX_GAP > "$LOG" 2>&1 && { cat "$LOG"; exit 0; }
  grep -qE "FileNotFoundError: .*site-packages|ModuleNotFoundError|ImportError" "$LOG" || { cat "$LOG"; exit 1; }
  sleep $((30 + RANDOM % 90))
done
cat "$LOG"; exit 1
