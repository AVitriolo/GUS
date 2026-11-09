#!/usr/bin/env bash
# Stability selection of the TCGA gene models, per cohort: a job array over the gene chunks (rebuild the train tables,
# selections on STAB_N_PAIRS complementary pairs of half-samples) + a summary job that starts when the array ends.
#   ./submit_stability.sh BRCA LUAD            chosen cohorts
#   ./submit_stability.sh all                  every cohort whose gene models are finished
#   ./submit_stability.sh --summary_only PRAD  only the summary (e.g. after changing STAB_PFER in the config)
#   ./submit_stability.sh --force PRAD         also cohorts whose gene models are not all finished
# Rerunning is safe: finished chunks are skipped. Progress: ./stability_status.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
if [ ! -f "$HERE/config/pipeline.conf" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then HERE=$SLURM_SUBMIT_DIR; fi
export PIPE_DIR=$HERE
source "$PIPE_DIR/config/pipeline.conf"
SUMMARY_ONLY=0; FORCE=0; COHORTS=()
for x in "$@"; do
  case $x in --summary_only) SUMMARY_ONLY=1 ;; --force) FORCE=1 ;; *) COHORTS+=("$x") ;; esac
done
[ ${#COHORTS[@]} -gt 0 ] || { echo "usage: $0 [--summary_only] [--force] COHORT [COHORT ...] | all"; exit 1; }
[ "${COHORTS[0]}" = all ] && COHORTS=($(cat "$WORK_DIR/cohorts_selected.txt" 2>/dev/null || ls "$WORK_DIR"))
EXTRA=""; [ -n "${SLURM_PARTITION:-}" ] && EXTRA="$EXTRA --partition=$SLURM_PARTITION"
[ -n "${SLURM_ACCOUNT:-}" ] && EXTRA="$EXTRA --account=$SLURM_ACCOUNT"; [ -n "${SLURM_QOS:-}" ] && EXTRA="$EXTRA --qos=$SLURM_QOS"

for COHORT in "${COHORTS[@]}"; do
  C=$WORK_DIR/$COHORT
  [ -f "$C/gene_list_$COHORT.txt" ] && ls $C/chunks/chunk_*.txt >/dev/null 2>&1 || { echo "$COHORT: no gene chunks - skipped"; continue; }
  tot=$(wc -l < "$C/gene_list_$COHORT.txt"); done=$(ls "$C/results_$METH_SCALE/performance" 2>/dev/null | wc -l)
  nod=$(cat "$C"/datasets_$METH_SCALE/*/status.tsv 2>/dev/null | awk '$2=="no_dataset"' | sort -u | wc -l)
  if [ $((done + nod)) -lt "$tot" ] && [ $FORCE = 0 ]; then
    echo "$COHORT: gene models not finished ($done + $nod without data of $tot) - skipped (--force to submit anyway)"; continue
  fi
  mkdir -p "$C/logs"
  DEP=""
  if [ $SUMMARY_ONLY = 0 ]; then
    N=$(ls $C/chunks/chunk_*.txt | wc -l)
    JID=$(sbatch --parsable --job-name=stab_$COHORT $EXTRA --array=0-$((N-1))%$STAB_MAX_CONCURRENT \
      --cpus-per-task=$STAB_CPUS --mem=$STAB_MEM --time=$STAB_TIME --output=$C/logs/stability_%A_%a.out \
      --export=ALL,PIPE_DIR=$PIPE_DIR,WORK_DIR=$WORK_DIR "$PIPE_DIR/slurm/stability.sbatch" $COHORT)
    DEP="--dependency=afterany:$JID"
    echo "$COHORT: stability array $JID ($N chunks)"
  fi
  SID=$(sbatch --parsable --job-name=stabsum_$COHORT $EXTRA $DEP --cpus-per-task=$STAB_SUM_CPUS --mem=$STAB_SUM_MEM \
    --time=$STAB_SUM_TIME --output=$C/logs/stability_summary_%j.out \
    --export=ALL,PIPE_DIR=$PIPE_DIR,WORK_DIR=$WORK_DIR "$PIPE_DIR/slurm/stability_summary.sbatch" $COHORT)
  echo "$COHORT: summary job $SID${DEP:+ (starts when the array ends)}"
done
