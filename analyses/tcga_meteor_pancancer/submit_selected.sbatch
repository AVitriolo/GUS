#!/usr/bin/env bash
#SBATCH --job-name=meteor_submit_selected
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/submit_selected_%j.out
# Submit the cohorts listed in $WORK_DIR/cohorts_selected.txt (written by download_local.sh). No internet needed:
# the prepare jobs unpack the copied .gz files. From the pipeline folder:  mkdir -p logs && sbatch submit_selected.sbatch
cd "$SLURM_SUBMIT_DIR"
export PIPE_DIR=$SLURM_SUBMIT_DIR; source config/pipeline.conf
[ -s "$WORK_DIR/cohorts_selected.txt" ] || { echo "no $WORK_DIR/cohorts_selected.txt - copy the downloaded folder first"; exit 1; }
echo "[$(date +%T)] submitting: $(tr '\n' ' ' < "$WORK_DIR/cohorts_selected.txt")"
for C in $(cat "$WORK_DIR/cohorts_selected.txt"); do
  ls "$WORK_DIR/$C/TCGA.$C.sampleMap/"HumanMethylation450* "$WORK_DIR/$C/TCGA.$C.sampleMap/"HiSeqV2* >/dev/null 2>&1 || echo "WARNING: $C data files missing in $WORK_DIR/$C"
done
bash ./submit_cohort.sh $(cat "$WORK_DIR/cohorts_selected.txt")
echo "[$(date +%T)] done - follow with ./status.sh"
