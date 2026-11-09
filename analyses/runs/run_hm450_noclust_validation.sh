#!/usr/bin/env bash
# Ablation of run_hm450_validation.sh: same HM450-probe-restricted CpGs, but skip
# merge_CpGs clustering entirely and use each surviving probe as its own feature
# ("simple probes"), instead of grouping them into metaCpG clusters ("CpGs in
# probes"). Tests the hypothesis that merge_CpGs_split_new.R's clustering grid
# (percentile bins 10-1000, tuned for TxIDs with hundreds of CpGs) is degenerate
# and actively harmful when only a handful of probe CpGs survive per TxID.
#
# Reuses data/rse_NTM_hm450/ and resources/CpG_data/top_K_hm450/ already produced
# by run_hm450_validation.sh for TxIDs that completed there (skips any TxID
# missing those, e.g. not yet processed by the main run) - does not touch or
# regenerate anything under data/rse_NTM_hm450, data/corr_NTM_hm450, or the
# clustered (_hm450, non-noclust) outputs.
#
# Usage: ./run_hm450_noclust_validation.sh TxID [TxID ...]
#   (unlike run_hm450_validation.sh, TxIDs are required here - this is meant to be
#   pointed at a specific set already completed by the main run, e.g. the first
#   200 "ok" TxIDs from logs/hm450_validation/_status_*.tsv)
#
# Parallelism / envs: same pattern as run_hm450_validation.sh - see that script's
# header comments for why PARALLEL_JOBS uses `xargs -P` and why *_ENV point at
# specific hashed .snakemake/conda/ directories rather than named conda envs.

set -uo pipefail

R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_  # matches workflow/envs/r_merge_CpGs.yml
PY_ML_ENV=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_         # matches workflow/envs/py_ML.yml

CV=3
N_ITER=27
N_JOBS=2
PARALLEL_JOBS="${PARALLEL_JOBS:-13}"

COUNTS_TRAIN=resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv
COUNTS_TEST=resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv
TRAIN_SAMPLES=resources/split/train_samples.txt
TEST_SAMPLES=resources/split/test_samples.txt
PREFIX=hg38_v38_pc_50000_25_10

mkdir -p data/xgb_noclust_hm450 \
         results/performance_noclust_hm450 results/sel_feats_noclust_hm450 \
         results/hyper_noclust_hm450 results/shap_noclust_hm450 \
         results/models logs/hm450_noclust_validation \
         resources/CpG_data/top_K_hm450_identity

STATUS_FILE="logs/hm450_noclust_validation/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local rse_dir="data/rse_NTM_hm450/${base}"
  local cpg_list="resources/CpG_data/top_K_hm450/${base}.bed"
  local log="logs/hm450_noclust_validation/${base}.log"

  if [ ! -d "$rse_dir" ] || [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  (
    echo "=== $TxID ==="

    # identity "mapping": each surviving probe CpG is its own cluster (no
    # clustering at all) - same 2-col CpGID<TAB>cluster_id format merge_CpGs
    # writes, so xgb_input_clustered_split.R needs no changes.
    local identity_map="resources/CpG_data/top_K_hm450_identity/${base}_identity_map.txt"
    awk -F'\t' '{print $2"\t"$2}' "$cpg_list" > "$identity_map"

    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="$rse_dir" \
      --input_path_mapping="$identity_map" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="data/xgb_noclust_hm450/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="data/xgb_noclust_hm450/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest.py \
      --input_path_dataset_train="data/xgb_noclust_hm450/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="data/xgb_noclust_hm450/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_noclust_hm450/${out_name}" \
      --output_path_performances="results/performance_noclust_hm450/${out_name}" \
      --output_path_selected_features="results/sel_feats_noclust_hm450/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_noclust_hm450/${out_name}"
  ) > "$log" 2>&1
  local status=$?

  if [ $status -ne 0 ]; then
    printf '%s\tfailed\n' "$TxID" >> "$STATUS_FILE"
  else
    printf '%s\tok\n' "$TxID" >> "$STATUS_FILE"
  fi
}
export -f process_txid run_r_merge_CpGs run_py_ML
export R_MERGE_CPGS_ENV PY_ML_ENV CV N_ITER N_JOBS \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES PREFIX STATUS_FILE

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 TxID [TxID ...]" >&2
  exit 1
fi
TxIDs=("$@")

echo "Running HM450 no-clustering (raw probes) validation for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_missing missing prerequisites (not yet done in main run), $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/hm450_noclust_validation/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
