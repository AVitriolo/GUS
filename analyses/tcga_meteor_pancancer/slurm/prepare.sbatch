#!/usr/bin/env bash
# Step 1 for one cohort: download (Xena), preprocess (REMP imputation, export), split the gene list into chunks,
# then submit the gene job array. Submitted by submit_cohort.sh (resources and log paths are set there).
set -euo pipefail
COHORT=$1
source "$PIPE_DIR/config/pipeline.conf"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
# BeeGFS sometimes fails to serve package files when many jobs start at once: start at a random moment, and retry an
# R step when it fails while LOADING packages ("cannot open file" / "namespace load failed"); real errors fail as usual
sleep $((RANDOM % 120))
retry_r() {   # retry_r <log> <command...>
  local log=$1; shift
  for try in 1 2 3; do
    rc=0; "$@" 2>&1 | tee "$log" || rc=$?      # "|| rc=" keeps set -e from ending the job before the retry check
    [ $rc -eq 0 ] && return 0
    grep -qE "cannot open file|namespace load failed|there is no package called" "$log" || return $rc
    echo "[$(date)] package loading failed (filesystem), retry $try/3 in 1-3 min"; sleep $((60 + RANDOM % 120))
  done
  return 1
}
C=$WORK_DIR/$COHORT; mkdir -p $C/chunks $C/logs
cd $C
# files downloaded elsewhere (download_local.sh) and copied here: unpack them; otherwise download (needs internet)
for f in HumanMethylation450 HiSeqV2; do
  g="TCGA.$COHORT.sampleMap/$f.gz"; u="TCGA.$COHORT.sampleMap/${f}_$COHORT"
  if [ ! -f "$u" ] && [ -f "$g" ]; then echo "[$(date)] $COHORT: unpacking $g"; gunzip -c "$g" > "$u.tmp" && mv "$u.tmp" "$u"; fi
done
echo "[$(date)] $COHORT: download"
[ -f "TCGA.$COHORT.sampleMap/HumanMethylation450_$COHORT" ] && [ -f "TCGA.$COHORT.sampleMap/HiSeqV2_$COHORT" ] \
  || retry_r $C/logs/download_R.log Rscript "$PIPE_DIR/scripts/01_download_xena.R" --cohort=$COHORT --out_dir=$C
echo "[$(date)] $COHORT: preprocess"
[ -f "need_gene_exp_$COHORT.csv" ] && [ -f "${COHORT}_data_meth_beta.csv" ] && [ -f "gene_list_$COHORT.txt" ] \
  || retry_r $C/logs/preprocess_R.log Rscript "$PIPE_DIR/scripts/02_preprocess.R" --cohort=$COHORT --work_dir=$C --resources="$PIPE_DIR/resources" --scripts="$PIPE_DIR/scripts" --cores=$SLURM_CPUS_PER_TASK
NB=$(head -1 "need_gene_exp_$COHORT.csv" | tr ',' '\n' | tail -n +2 | wc -l)
echo "[$(date)] $COHORT: $NB samples with both methylation and expression after preprocessing"
if [ "$NB" -lt "$MIN_SAMPLES" ]; then echo "$COHORT: fewer than $MIN_SAMPLES matched samples - not modelled"; exit 0; fi
echo "[$(date)] $COHORT: chunks"
rm -f chunks/chunk_*.txt; split -l $GENES_PER_TASK -d -a 4 --additional-suffix=.txt "gene_list_$COHORT.txt" chunks/chunk_
N=$(ls chunks/chunk_*.txt | wc -l); echo "$(wc -l < gene_list_$COHORT.txt) genes in $N chunks of $GENES_PER_TASK"
EXTRA=""; [ -n "${SLURM_PARTITION:-}" ] && EXTRA="$EXTRA --partition=$SLURM_PARTITION"
[ -n "${SLURM_ACCOUNT:-}" ] && EXTRA="$EXTRA --account=$SLURM_ACCOUNT"; [ -n "${SLURM_QOS:-}" ] && EXTRA="$EXTRA --qos=$SLURM_QOS"
sbatch --job-name=meteor_$COHORT $EXTRA --array=0-$((N-1))%$MAX_CONCURRENT_TASKS \
  --cpus-per-task=$((PAR_PER_TASK * N_JOBS)) --mem=$TASK_MEM --time=$TASK_TIME \
  --output=$C/logs/genes_%A_%a.out --export=ALL,PIPE_DIR=$PIPE_DIR,WORK_DIR=$WORK_DIR "$PIPE_DIR/slurm/genes.sbatch" $COHORT
echo "[$(date)] $COHORT: gene array submitted"
