#!/usr/bin/env bash
# Same question as run_cv_stability_validation.sh (does RepeatedKFold fix cv_r2's
# instability?), but on the UNRESTRICTED CpG set: every WGBS CpG within the 2Mb
# window (+-1Mb of TSS), not just the ones that are also HM450 probes. Isolates
# the effect of probe-restriction alone, at a window size where that's no longer
# confounded by window width (both this and the HM450 2Mb condition use the
# identical +-1Mb window).
#
# Candidate lists (resources/CpG_data/top_K_2Mb_full/{base}.bed) were built from
# the SAME wgbs_cpgs_in_window intermediate used for the HM450 2Mb candidates,
# just without the probe-intersection step - ~14,681 CpGs/TxID on average vs
# ~175 for the probe-restricted version.
#
# No clustering (matches the HM450 2Mb no-clustering condition this compares
# against) - each CpG stays its own feature. Uses model_selection_split_kbest_repeated.py
# (RepeatedKFold) directly, since there's no non-repeated baseline worth running
# first here - cv=3/n_repeats=5 already established as the useful setting.
#
# Usage: ./run_2mb_full_cv_stability_validation.sh TxID [TxID ...]
#
# THREAD SAFETY: a first attempt at this (unrestricted, ~14,681 CpGs/TxID) blew
# the machine's load to 249 on 112 cores at PARALLEL_JOBS=40, because R's BLAS
# and sklearn/joblib's worker pool both multi-thread internally BY DEFAULT,
# uncapped - so real concurrency was PARALLEL_JOBS * N_JOBS * (however many BLAS
# threads each process grabbed), not PARALLEL_JOBS * N_JOBS. The env vars below
# pin every BLAS/OpenMP backend to 1 thread per process, so PARALLEL_JOBS *
# N_JOBS is actually what you get. Keep PARALLEL_JOBS in the 10-15 range even
# with this fix - this condition's R/SVR fits are individually much heavier
# (~14,681 vs ~175 features/TxID) than anything else run in this investigation.

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
         results/performance_2mb_full_cv3r5 results/sel_feats_2mb_full_cv3r5 \
         results/hyper_2mb_full_cv3r5 results/shap_2mb_full_cv3r5 \
         results/models logs/2mb_full_cv_stability \
         resources/CpG_data/top_K_2Mb_full_identity

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 TxID [TxID ...]" >&2
  exit 1
fi
TxIDs=("$@")

STATUS_FILE="logs/2mb_full_cv_stability/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb_full/${base}.bed"
  local log="logs/2mb_full_cv_stability/${base}.log"

  if [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="logs/2mb_full_cv_stability/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ($(wc -l < "$cpg_list") candidate CpGs) ==="

    run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
      --input_dir="$INTEGRATED_H5" \
      --input_path_filtered_CpGs="$cpg_list" \
      --output_path_rse="data/rse_NTM_2mb_full/${base}" \
    || exit 1

    local identity_map="resources/CpG_data/top_K_2Mb_full_identity/${base}_identity_map.txt"
    awk -F'\t' '{print $2"\t"$2}' "$cpg_list" > "$identity_map"

    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="data/rse_NTM_2mb_full/${base}" \
      --input_path_mapping="$identity_map" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="data/xgb_2mb_full/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="data/xgb_2mb_full/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest_repeated.py \
      --input_path_dataset_train="data/xgb_2mb_full/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="data/xgb_2mb_full/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --n_repeats="$N_REPEATS" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_2mb_full_cv3r5/${out_name}" \
      --output_path_performances="results/performance_2mb_full_cv3r5/${out_name}" \
      --output_path_selected_features="results/sel_feats_2mb_full_cv3r5/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_2mb_full_cv3r5/${out_name}" \
      --skip_shap
  ) > "$log" 2>&1
  local status=$?
  rm -rf "$txid_tmp"

  # Clean up only THIS TxID's own large intermediates now that
  # model_selection is done with them - does not touch any other TxID's
  # files and does not skip/cancel any processing. At full 27,648-TxID
  # scale these are not safe to accumulate unbounded (this run filled a
  # 3.6TB shared filesystem to 100% partway through: ~34MB RSE + a wide
  # feature-matrix pair + a pickled model per TxID, times thousands of
  # TxIDs). Only the small performance/sel_feats/hyper text outputs are kept.
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

echo "Running 2Mb FULL-CpG CV-stability validation (cv=$CV, n_repeats=$N_REPEATS) for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_missing missing prerequisites, $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/2mb_full_cv_stability/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
