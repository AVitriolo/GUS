#!/usr/bin/env bash
#SBATCH --job-name=meteor_genelists
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/rebuild_gene_lists_%j.out
# For every selected cohort whose preparation is finished: rebuild the gene list with all expressed genes.
# From the pipeline folder:  sbatch slurm/rebuild_gene_lists.sbatch   (then: sbatch submit_selected.sbatch)
cd "$SLURM_SUBMIT_DIR"; export PIPE_DIR=$SLURM_SUBMIT_DIR; source config/pipeline.conf
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
for C in $(cat "$WORK_DIR/cohorts_selected.txt"); do
  if [ -f "$WORK_DIR/$C/$C.data.Rdata" ] && [ -f "$WORK_DIR/$C/gene_list_$C.txt" ]; then
    Rscript scripts/rebuild_gene_list.R --cohort=$C --work_dir="$WORK_DIR/$C"
  else echo "$C: preparation not finished yet - rerun this job later"; fi
done
