#!/usr/bin/env bash
# Step 2: one array task = one chunk of genes (chunks/chunk_<task id>.txt). Submitted by prepare.sbatch.
set -uo pipefail
COHORT=$1
source "$PIPE_DIR/config/pipeline.conf"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
sleep $((RANDOM % 120))   # spread the start of the array tasks (BeeGFS)
CHUNK=$(printf "%s/%s/chunks/chunk_%04d.txt" "$WORK_DIR" "$COHORT" "$SLURM_ARRAY_TASK_ID")
[ -f "$CHUNK" ] || { echo "no chunk $CHUNK"; exit 1; }
bash "$PIPE_DIR/scripts/run_gene_chunk.sh" $COHORT "$CHUNK"
