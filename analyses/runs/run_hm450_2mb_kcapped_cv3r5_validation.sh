#!/usr/bin/env bash
# Same candidate-CpG pipeline as run_hm450_2mb_validation.sh (HM450-probe-
# restricted, 2Mb TSS window, "simple probes as individual features, no
# clustering"), but with the modeling step swapped for
# model_selection_split_kbest_repeated_kcapped.py (RepeatedKFold cv=3/
# n_repeats=5, k capped at per-fold training size) - the same modeling script
# and parameters used by run_2mb_full_kcapped_cv_stability_validation.sh for
# the unrestricted-CpG "WGBS complete" baseline, so the two are directly
# comparable (same window, same modeling protocol, only the candidate CpG set
# differs).
#
# Requires resources/CpG_data/top_K_2Mb/{base}.bed to already exist for each
# TxID (real WGBS CpGs in the 2Mb window, filtered down to HM450 probes - see
# the saved memory cpg-candidate-set-construction-order.md for why that order
# matters). This script does not regenerate them.
#
# Usage: ./run_hm450_2mb_kcapped_cv3r5_validation.sh TxID [TxID ...]

set -uo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

R_XGB_INPUT_ENV=/home/Prostate_GAS/.snakemake/conda/1be302004981d9afc5efa19c4056445a_   # matches workflow/envs/r_xgb_input.yml
R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_  # matches workflow/envs/r_merge_CpGs.yml
PY_ML_ENV=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_         # matches workflow/envs/py_ML.yml

CV=3
N_REPEATS=5
N_ITER=27
N_JOBS=2
PARALLEL_JOBS="${PARALLEL_JOBS:-12}"

COUNTS_TRAIN=resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv
COUNTS_TEST=resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv
TRAIN_SAMPLES=resources/split/train_samples.txt
TEST_SAMPLES=resources/split/test_samples.txt
INTEGRATED_H5=resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5
PREFIX=hg38_v38_pc_50000_25_10

mkdir -p data/rse_NTM_2mb_hm450_kcapped data/xgb_2mb_hm450_kcapped \
         results/performance_2mb_hm450_kcapped_cv3r5 results/sel_feats_2mb_hm450_kcapped_cv3r5 \
         results/hyper_2mb_hm450_kcapped_cv3r5 results/shap_2mb_hm450_kcapped_cv3r5 \
         results/models logs/2mb_hm450_kcapped_cv_stability \
         resources/CpG_data/top_K_2Mb_identity

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 TxID [TxID ...]" >&2
  exit 1
fi
TxIDs=("$@")

STATUS_FILE="logs/2mb_hm450_kcapped_cv_stability/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb/${base}.bed"
  local perf_file="results/performance_2mb_hm450_kcapped_cv3r5/${PREFIX}_${CV}_${TxID}"
  local log="logs/2mb_hm450_kcapped_cv_stability/${base}.log"

  if [ -f "$perf_file" ]; then
    printf '%s\talready_done\n' "$TxID" >> "$STATUS_FILE"
    return
  fi
  if [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="logs/2mb_hm450_kcapped_cv_stability/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ==="

    # 1. RSE for the HM450-restricted 2Mb candidate probe set
    run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
      --input_dir="$INTEGRATED_H5" \
      --input_path_filtered_CpGs="$cpg_list" \
      --output_path_rse="data/rse_NTM_2mb_hm450_kcapped/${base}" \
    || exit 1

    # 2. identity "mapping" - each probe is its own feature, no clustering
    local identity_map="resources/CpG_data/top_K_2Mb_identity/${base}_identity_map.txt"
    awk -F'\t' '{print $2"\t"$2}' "$cpg_list" > "$identity_map"

    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="data/rse_NTM_2mb_hm450_kcapped/${base}" \
      --input_path_mapping="$identity_map" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="data/xgb_2mb_hm450_kcapped/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="data/xgb_2mb_hm450_kcapped/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    # 3. model - RepeatedKFold cv=3/n_repeats=5, k capped (matches the WGBS
    # complete kcapped baseline's modeling protocol exactly)
    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML workflow/scripts/model_selection_split_kbest_repeated_kcapped.py \
      --input_path_dataset_train="data/xgb_2mb_hm450_kcapped/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="data/xgb_2mb_hm450_kcapped/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --n_repeats="$N_REPEATS" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_2mb_hm450_kcapped_cv3r5/${out_name}" \
      --output_path_performances="$perf_file" \
      --output_path_selected_features="results/sel_feats_2mb_hm450_kcapped_cv3r5/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_2mb_hm450_kcapped_cv3r5/${out_name}" \
      --skip_shap
  ) > "$log" 2>&1
  local status=$?

  # Clean up large per-TxID intermediates - not safe to accumulate unbounded
  # at ~23k-TxID scale (same rationale as run_2mb_full_kcapped_cv_stability_validation.sh).
  rm -rf "data/rse_NTM_2mb_hm450_kcapped/${base}"
  rm -f "data/xgb_2mb_hm450_kcapped/${base}_xgb_input_clustered_train.txt" \
        "data/xgb_2mb_hm450_kcapped/${base}_xgb_input_clustered_test.txt"
  rm -f results/models/${TxID}_*.pkl
  rm -rf "$txid_tmp"

  if [ $status -ne 0 ]; then
    printf '%s\tfailed\n' "$TxID" >> "$STATUS_FILE"
  else
    printf '%s\tok\n' "$TxID" >> "$STATUS_FILE"
  fi
}
export -f process_txid run_r_xgb_input run_r_merge_CpGs run_py_ML
export R_XGB_INPUT_ENV R_MERGE_CPGS_ENV PY_ML_ENV CV N_REPEATS N_ITER N_JOBS \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES INTEGRATED_H5 PREFIX STATUS_FILE

echo "Running HM450 2Mb-window k-CAPPED CV-stability validation (cv=$CV, n_repeats=$N_REPEATS) for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_done=$(awk -F'\t' '$2=="already_done"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_done already done, $n_missing missing prerequisites, $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/2mb_hm450_kcapped_cv_stability/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
