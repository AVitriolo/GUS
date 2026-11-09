#!/usr/bin/env bash
# Companion to run_hm450_2mb_validation.sh: same 2Mb-window candidate probe pool
# (matching TCGA PRAD's actual +-1,000,000bp window), but WITH merge_CpGs
# clustering (unlike the no-clustering ablations). Tests whether clustering, which
# was degenerate/harmful at the tiny 50kb-window candidate counts (10-20 CpGs),
# behaves reasonably once there's a much richer candidate pool per TxID (~150-200
# candidates at 2Mb) - i.e. whether "simple probes" beats "CpGs in probes" only
# because of scarcity, or whether no-clustering is just better outright even with
# plenty of CpGs to cluster.
#
# Uses its own RSE directory (data/rse_NTM_2mb_clustered/, not
# data/rse_NTM_2mb/) even though the candidate CpG list is identical to
# run_hm450_2mb_validation.sh's - deliberately avoids sharing that directory so
# this can run concurrently with the still-running no-clustering 2Mb job without
# two processes writing the same HDF5Array output directory at once.
#
# For each TxID: unpack_RSE_NTM.R -> compute_pairwise_correlation_split.R (new,
# needed - no existing 2Mb correlation matrix exists anywhere) -> merge_CpGs_split_new.R
# -> xgb_input_clustered_split.R -> model_selection_split_kbest.py. All existing,
# unmodified scripts.
#
# Usage: ./run_hm450_2mb_clustered_validation.sh TxID [TxID ...]

set -uo pipefail

R_XGB_INPUT_ENV=/home/Prostate_GAS/.snakemake/conda/1be302004981d9afc5efa19c4056445a_   # matches workflow/envs/r_xgb_input.yml
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
INTEGRATED_H5=resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5
PREFIX=hg38_v38_pc_50000_25_10

mkdir -p data/rse_NTM_2mb_clustered data/corr_NTM_2mb \
         data/clustering/clusters_maps_2mb results/plots/clustering_2mb \
         data/xgb_2mb_clustered \
         results/performance_2mb_clustered results/sel_feats_2mb_clustered \
         results/hyper_2mb_clustered results/shap_2mb_clustered \
         results/models logs/hm450_2mb_clustered_validation

STATUS_FILE="logs/hm450_2mb_clustered_validation/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb/${base}.bed"
  local log="logs/hm450_2mb_clustered_validation/${base}.log"

  if [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="logs/hm450_2mb_clustered_validation/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ==="

    # 1. RSE for the 2Mb candidate probe set (own directory, see header note)
    run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
      --input_dir="$INTEGRATED_H5" \
      --input_path_filtered_CpGs="$cpg_list" \
      --output_path_rse="data/rse_NTM_2mb_clustered/${base}" \
    || exit 1

    # 2. pairwise correlation matrix (needed for clustering - doesn't exist yet at 2Mb)
    run_r_merge_CpGs workflow/scripts/compute_pairwise_correlation_split.R \
      --input_path_rse="data/rse_NTM_2mb_clustered/${base}" \
      --train_samples="$TRAIN_SAMPLES" \
      --output_path_corr="data/corr_NTM_2mb/${base}" \
      --tx_id="$TxID" \
    || exit 1

    # 3. cluster
    run_r_merge_CpGs workflow/scripts/merge_CpGs_split_new.R \
      --input_path="data/corr_NTM_2mb/${base}" \
      --input_path_rse="data/rse_NTM_2mb_clustered/${base}" \
      --input_path_counts="$COUNTS_TRAIN" \
      --train_samples="$TRAIN_SAMPLES" \
      --output_path_plot="results/plots/clustering_2mb/${base}_heatmap.pdf" \
      --output_path_map="data/clustering/clusters_maps_2mb/${base}_clusters_map.txt" \
      --C=2 --j=1 \
    || exit 1

    # 4. clustered feature matrix
    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="data/rse_NTM_2mb_clustered/${base}" \
      --input_path_mapping="data/clustering/clusters_maps_2mb/${base}_clusters_map.txt" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="data/xgb_2mb_clustered/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="data/xgb_2mb_clustered/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    # 5. model
    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest.py \
      --input_path_dataset_train="data/xgb_2mb_clustered/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="data/xgb_2mb_clustered/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_2mb_clustered/${out_name}" \
      --output_path_performances="results/performance_2mb_clustered/${out_name}" \
      --output_path_selected_features="results/sel_feats_2mb_clustered/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_2mb_clustered/${out_name}"
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
export R_XGB_INPUT_ENV R_MERGE_CPGS_ENV PY_ML_ENV CV N_ITER N_JOBS \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES INTEGRATED_H5 PREFIX STATUS_FILE

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 TxID [TxID ...]" >&2
  exit 1
fi
TxIDs=("$@")

echo "Running HM450 2Mb-window CLUSTERED validation for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_missing missing prerequisites, $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/hm450_2mb_clustered_validation/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
