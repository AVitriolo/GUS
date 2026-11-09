#!/usr/bin/env bash
# Lambda-sweep ablation of run_hm450_2mb_clustered_validation.sh: same 2Mb-window
# candidate probe pool and same clustering pipeline, but merge_CpGs_split_new_fragpen.R
# in place of merge_CpGs_split_new.R (additive, independently-weighted fragmentation
# penalty - see that script's header comment for why). Tests whether this reduces
# over-fragmented clustering and improves cv_r2/test_r2 vs the original formula.
#
# Builds its own RSE (data/rse_NTM_2mb_fragpen/) + correlation matrix
# (data/corr_NTM_2mb_fragpen/) cache, SHARED across lambda values (skip-if-exists,
# since RSE/corr don't depend on lambda - only the clustering step does), so a
# second/third lambda run only re-pays for clustering + xgb_input + model_selection,
# not the expensive h5 read + correlation computation.
#
# Usage: ./run_hm450_2mb_fragpen_validation.sh LAMBDA TxID [TxID ...]

set -uo pipefail

R_XGB_INPUT_ENV=/home/Prostate_GAS/.snakemake/conda/1be302004981d9afc5efa19c4056445a_   # matches workflow/envs/r_xgb_input.yml
R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_  # matches workflow/envs/r_merge_CpGs.yml
PY_ML_ENV=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_         # matches workflow/envs/py_ML.yml

CV=3
N_ITER=27
N_JOBS=2
PARALLEL_JOBS="${PARALLEL_JOBS:-8}"

COUNTS_TRAIN=resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv
COUNTS_TEST=resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv
TRAIN_SAMPLES=resources/split/train_samples.txt
TEST_SAMPLES=resources/split/test_samples.txt
INTEGRATED_H5=resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5
PREFIX=hg38_v38_pc_50000_25_10

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 LAMBDA TxID [TxID ...]" >&2
  exit 1
fi
LAMBDA="$1"; shift
TxIDs=("$@")
LTAG="L${LAMBDA}"

CLUST_DIR="data/clustering/clusters_maps_2mb_fragpen_${LTAG}"
PLOT_DIR="results/plots/clustering_2mb_fragpen_${LTAG}"
XGB_DIR="data/xgb_2mb_fragpen_${LTAG}"
PERF_DIR="results/performance_2mb_fragpen_${LTAG}"
SELF_DIR="results/sel_feats_2mb_fragpen_${LTAG}"
HYPER_DIR="results/hyper_2mb_fragpen_${LTAG}"
SHAP_DIR="results/shap_2mb_fragpen_${LTAG}"
LOG_DIR="logs/hm450_2mb_fragpen_validation_${LTAG}"

mkdir -p data/rse_NTM_2mb_fragpen data/corr_NTM_2mb_fragpen \
         "$CLUST_DIR" "$PLOT_DIR" "$XGB_DIR" \
         "$PERF_DIR" "$SELF_DIR" "$HYPER_DIR" "$SHAP_DIR" \
         results/models "$LOG_DIR"

STATUS_FILE="${LOG_DIR}/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb/${base}.bed"
  local log="${LOG_DIR}/${base}.log"

  if [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="${LOG_DIR}/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID (lambda=$LAMBDA) ==="

    local rse_dir="data/rse_NTM_2mb_fragpen/${base}"
    local corr_path="data/corr_NTM_2mb_fragpen/${base}"

    # 1. RSE - shared across lambda values (and across concurrently-running lambda
    # sweeps). mkdir is atomic, so use a lock dir to make sure only one concurrent
    # process builds a given TxID's cache; others wait for it instead of racing.
    if [ ! -d "$rse_dir" ]; then
      local rse_lock="${rse_dir}.building"
      if mkdir "$rse_lock" 2>/dev/null; then
        run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
          --input_dir="$INTEGRATED_H5" \
          --input_path_filtered_CpGs="$cpg_list" \
          --output_path_rse="$rse_dir" \
        || { rmdir "$rse_lock"; exit 1; }
        rmdir "$rse_lock"
      else
        echo "[$TxID] waiting for concurrent RSE build..."
        while [ -d "$rse_lock" ]; do sleep 2; done
        [ -d "$rse_dir" ] || { echo "[$TxID] concurrent RSE build did not succeed"; exit 1; }
      fi
    else
      echo "[$TxID] reusing existing RSE"
    fi

    # 2. correlation matrix - same shared-cache locking pattern
    if [ ! -f "$corr_path" ]; then
      local corr_lock="${corr_path}.building"
      if mkdir "$corr_lock" 2>/dev/null; then
        run_r_merge_CpGs workflow/scripts/compute_pairwise_correlation_split.R \
          --input_path_rse="$rse_dir" \
          --train_samples="$TRAIN_SAMPLES" \
          --output_path_corr="$corr_path" \
          --tx_id="$TxID" \
        || { rmdir "$corr_lock"; exit 1; }
        rmdir "$corr_lock"
      else
        echo "[$TxID] waiting for concurrent correlation build..."
        while [ -d "$corr_lock" ]; do sleep 2; done
        [ -f "$corr_path" ] || { echo "[$TxID] concurrent correlation build did not succeed"; exit 1; }
      fi
    else
      echo "[$TxID] reusing existing correlation matrix"
    fi

    # 3. cluster with the fragmentation-penalty variant
    run_r_merge_CpGs analyses/scripts/merge_CpGs_split_new_fragpen.R \
      --input_path="$corr_path" \
      --input_path_rse="$rse_dir" \
      --input_path_counts="$COUNTS_TRAIN" \
      --train_samples="$TRAIN_SAMPLES" \
      --output_path_plot="${PLOT_DIR}/${base}_heatmap.pdf" \
      --output_path_map="${CLUST_DIR}/${base}_clusters_map.txt" \
      --C=2 --j=1 --lambda="$LAMBDA" \
    || exit 1

    # 4. clustered feature matrix
    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="$rse_dir" \
      --input_path_mapping="${CLUST_DIR}/${base}_clusters_map.txt" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="${XGB_DIR}/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="${XGB_DIR}/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    # 5. model
    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest.py \
      --input_path_dataset_train="${XGB_DIR}/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="${XGB_DIR}/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="${HYPER_DIR}/${out_name}" \
      --output_path_performances="${PERF_DIR}/${out_name}" \
      --output_path_selected_features="${SELF_DIR}/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="${SHAP_DIR}/${out_name}"
  ) > "$log" 2>&1
  local status=$?
  rm -rf "$txid_tmp"

  if [ $status -ne 0 ]; then
    printf '%s\tfailed\n' "$TxID" >> "$STATUS_FILE"
  else
    printf '%s\tok\n' "$TxID" >> "$STATUS_FILE"
  fi
}
export -f process_txid run_r_xgb_input run_r_merge_CpGs run_py_ML
export R_XGB_INPUT_ENV R_MERGE_CPGS_ENV PY_ML_ENV CV N_ITER N_JOBS LAMBDA \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES INTEGRATED_H5 PREFIX \
       STATUS_FILE CLUST_DIR PLOT_DIR XGB_DIR PERF_DIR SELF_DIR HYPER_DIR SHAP_DIR LOG_DIR

echo "Running HM450 2Mb fragmentation-penalty (lambda=$LAMBDA) validation for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_missing missing prerequisites, $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see ${LOG_DIR}/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
