#!/usr/bin/env bash
# Stability summary for one cohort, after its stability array (submit_stability.sh sets the dependency):
#   stability_selection_summary.R  per feature: selection frequency, stable (PFER <= STAB_PFER), Nogueira phi per gene
#   stability_selection_regions.R  per co-methylated region: frequency, stable, phi per gene
#   stability_report.py            per-gene table (with cv / test R2) and plot by CV R2 bin
# To change STAB_PFER later, edit the config and resubmit only this job:  ./submit_stability.sh --summary_only COHORT
set -euo pipefail
COHORT=$1
source "$PIPE_DIR/config/pipeline.conf"
source "$(conda info --base)/etc/profile.d/conda.sh"
C=$WORK_DIR/$COHORT; S=$C/stability_$METH_SCALE
NC=$(ls $C/chunks/chunk_*.txt | wc -l); ND=$(ls $S/done_chunk_* 2>/dev/null | wc -l)
echo "[$(date)] $COHORT: $ND of $NC chunks done"
[ "$ND" -lt "$NC" ] && echo "WARNING: $((NC - ND)) chunks missing - their genes are left out (rerun ./submit_stability.sh $COHORT)"
conda activate "$STAB_R_ENV"
echo "[$(date)] $COHORT: feature-level stability"
Rscript "$PIPE_DIR/scripts/stability_selection_summary.R" "$S" "$STAB_N_PAIRS" "$STAB_PFER"
echo "[$(date)] $COHORT: region-level stability"
Rscript "$PIPE_DIR/scripts/stability_selection_regions.R" "$S" none "$STAB_MAX_GAP" "$STAB_N_PAIRS" "$STAB_PFER" \
  "$SLURM_CPUS_PER_TASK" "corr$STAB_CORR"
conda activate "$CONDA_ENV"
echo "[$(date)] $COHORT: report"
python "$PIPE_DIR/scripts/stability_report.py" --stab_dir "$S" --perf_dir "$C/results_$METH_SCALE/performance" \
  --tag "corr$STAB_CORR" --label "TCGA $COHORT"
echo "[$(date)] $COHORT: done -> $S/stability_report.png, stability_report_genes.tsv, cpg_freq.tsv, region_freq_corr$STAB_CORR.tsv"
