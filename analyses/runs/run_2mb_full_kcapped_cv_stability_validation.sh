#!/usr/bin/env bash
# Same as run_2mb_full_cv_stability_validation.sh (unrestricted 2Mb-window CpG
# set, RepeatedKFold cv=3/n_repeats=5), except uses
# model_selection_split_kbest_repeated_kcapped.py, which caps
# feature_selection__k at the per-fold training size instead of leaving it
# uncapped up to n_feat (~14,681). The uncapped run picked k=11,540 for one
# TxID - see GET_HYPERPARAM_SPACE_KBEST_capped.py's header for the full
# reasoning.
#
# Shares the RSE (data/rse_NTM_2mb_full/) and xgb feature-matrix
# (data/xgb_2mb_full/) caches with the uncapped run via mkdir-lock guarded
# reuse (same pattern as run_hm450_2mb_fragpen_validation.sh's lambda sweep) -
# neither of those steps depends on the k-cap, only the modeling step does, so
# this can run SAFELY CONCURRENTLY with the still-in-flight uncapped batch
# without duplicate builds or write races.
#
# Usage: ./run_2mb_full_kcapped_cv_stability_validation.sh TxID [TxID ...]

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

mkdir -p data/rse_NTM_2mb_full data/xgb_2mb_full \
         results/performance_2mb_full_kcapped_cv3r5 results/sel_feats_2mb_full_kcapped_cv3r5 \
         results/hyper_2mb_full_kcapped_cv3r5 results/shap_2mb_full_kcapped_cv3r5 \
         results/models logs/2mb_full_kcapped_cv_stability \
         resources/CpG_data/top_K_2Mb_full_identity

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 TxID [TxID ...]" >&2
  exit 1
fi
TxIDs=("$@")

STATUS_FILE="logs/2mb_full_kcapped_cv_stability/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb_full/${base}.bed"
  local log="logs/2mb_full_kcapped_cv_stability/${base}.log"

  if [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="logs/2mb_full_kcapped_cv_stability/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ==="

    local rse_dir="data/rse_NTM_2mb_full/${base}"

    # 1. RSE - shared with the uncapped run, mkdir-lock guarded
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

    # 2. clustered feature matrix (identity mapping, no clustering) - shared, same lock pattern
    local train_xgb="data/xgb_2mb_full/${base}_xgb_input_clustered_train.txt"
    local test_xgb="data/xgb_2mb_full/${base}_xgb_input_clustered_test.txt"
    if [ ! -f "$train_xgb" ] || [ ! -f "$test_xgb" ]; then
      local xgb_lock="data/xgb_2mb_full/${base}.building"
      if mkdir "$xgb_lock" 2>/dev/null; then
        local identity_map="resources/CpG_data/top_K_2Mb_full_identity/${base}_identity_map.txt"
        awk -F'\t' '{print $2"\t"$2}' "$cpg_list" > "$identity_map"

        run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
          --input_path_rse="$rse_dir" \
          --input_path_mapping="$identity_map" \
          --input_path_counts_train="$COUNTS_TRAIN" \
          --input_path_counts_test="$COUNTS_TEST" \
          --train_samples="$TRAIN_SAMPLES" \
          --test_samples="$TEST_SAMPLES" \
          --output_path_train="$train_xgb" \
          --output_path_test="$test_xgb" \
        || { rmdir "$xgb_lock"; exit 1; }
        rmdir "$xgb_lock"
      else
        echo "[$TxID] waiting for concurrent xgb_input build..."
        while [ -d "$xgb_lock" ]; do sleep 2; done
        [ -f "$train_xgb" ] && [ -f "$test_xgb" ] || { echo "[$TxID] concurrent xgb_input build did not succeed"; exit 1; }
      fi
    else
      echo "[$TxID] reusing existing xgb feature matrix"
    fi

    # 3. model - the only step that actually differs from the uncapped run
    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML workflow/scripts/model_selection_split_kbest_repeated_kcapped.py \
      --input_path_dataset_train="$train_xgb" \
      --input_path_dataset_test="$test_xgb" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --n_repeats="$N_REPEATS" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_2mb_full_kcapped_cv3r5/${out_name}" \
      --output_path_performances="results/performance_2mb_full_kcapped_cv3r5/${out_name}" \
      --output_path_selected_features="results/sel_feats_2mb_full_kcapped_cv3r5/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_2mb_full_kcapped_cv3r5/${out_name}" \
      --skip_shap
  ) > "$log" 2>&1
  local status=$?
  rm -rf "$txid_tmp"

  # Clean up this TxID's large intermediates now that model_selection is done
  # with them. NOTE: the "shared cache" reuse above (mkdir-lock guarded) is
  # now mostly theoretical for concurrent runs against the uncapped script,
  # since that script ALSO cleans up its own RSE/xgb immediately after use
  # (added after a prior run filled a 3.6TB shared filesystem to 100%) - so
  # there's usually nothing left to reuse in practice. Clean up here too, for
  # the same reason: these are not safe to accumulate unbounded at
  # 27,648-TxID scale. Only the small performance/sel_feats/hyper text
  # outputs are kept.
  rm -rf "data/rse_NTM_2mb_full/${base}"
  rm -f "data/xgb_2mb_full/${base}_xgb_input_clustered_train.txt" \
        "data/xgb_2mb_full/${base}_xgb_input_clustered_test.txt"
  rm -f results/models/${TxID}_*.pkl

  if [ $status -ne 0 ]; then
    printf '%s\tfailed\n' "$TxID" >> "$STATUS_FILE"
  else
    printf '%s\tok\n' "$TxID" >> "$STATUS_FILE"
  fi
}
export -f process_txid run_r_xgb_input run_r_merge_CpGs run_py_ML
export R_XGB_INPUT_ENV R_MERGE_CPGS_ENV PY_ML_ENV CV N_REPEATS N_ITER N_JOBS \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES INTEGRATED_H5 PREFIX STATUS_FILE

echo "Running 2Mb FULL-CpG k-CAPPED CV-stability validation (cv=$CV, n_repeats=$N_REPEATS) for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_missing missing prerequisites, $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/2mb_full_kcapped_cv_stability/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
