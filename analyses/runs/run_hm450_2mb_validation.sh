#!/usr/bin/env bash
# Third ablation, widening the candidate-probe window further to match TCGA PRAD's
# actual methodology exactly: /home/ieo7293/TCGA_PRAD/launch_gus_prad.sh calls
# make_datasets_prad_positional_cli.py --window_bp 1000000, i.e. +-1,000,000bp from
# the TSS-anchor probe = a 2Mb TOTAL window - twice as wide as the 1Mb-total window
# tested in run_hm450_1mb_validation.sh. Same "simple probes as individual features,
# no clustering" approach as that script.
#
# Candidate probe lists (resources/CpG_data/top_K_2Mb/{base}.bed, 2-col
# TxID<TAB>CpGID) must already exist - generated once via bedtools slop(-b 1000000)
# + intersect, in the CORRECT order (real WGBS CpGs in the window first, filtered
# down to Illumina probes second - see the saved memory
# cpg-candidate-set-construction-order.md and the plan file for the exact commands
# used). This script does not regenerate them.
#
# For each TxID: regenerates the RSE from scratch via the existing, unmodified
# unpack_RSE_NTM.R, builds an identity mapping (each probe its own feature), then
# the existing, unmodified xgb_input_clustered_split.R and model_selection_split_kbest.py.
#
# Note: PREFIX stays "hg38_v38_pc_50000_25_10" for filename consistency with the
# other conditions, even though the candidate window here is 2Mb, not the 50kb
# baked into that historical string.
#
# Usage: ./run_hm450_2mb_validation.sh TxID [TxID ...]

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

mkdir -p data/rse_NTM_2mb data/xgb_2mb_hm450 \
         results/performance_2mb_hm450 results/sel_feats_2mb_hm450 \
         results/hyper_2mb_hm450 results/shap_2mb_hm450 \
         results/models logs/hm450_2mb_validation \
         resources/CpG_data/top_K_2Mb_identity

STATUS_FILE="logs/hm450_2mb_validation/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb/${base}.bed"
  local log="logs/hm450_2mb_validation/${base}.log"

  if [ ! -f "$cpg_list" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="logs/hm450_2mb_validation/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ==="

    # 1. regenerate the RSE for the 2Mb candidate probe set
    run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
      --input_dir="$INTEGRATED_H5" \
      --input_path_filtered_CpGs="$cpg_list" \
      --output_path_rse="data/rse_NTM_2mb/${base}" \
    || exit 1

    # 2. identity "mapping" - each probe is its own feature, no clustering
    local identity_map="resources/CpG_data/top_K_2Mb_identity/${base}_identity_map.txt"
    awk -F'\t' '{print $2"\t"$2}' "$cpg_list" > "$identity_map"

    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="data/rse_NTM_2mb/${base}" \
      --input_path_mapping="$identity_map" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="data/xgb_2mb_hm450/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="data/xgb_2mb_hm450/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest.py \
      --input_path_dataset_train="data/xgb_2mb_hm450/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="data/xgb_2mb_hm450/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_2mb_hm450/${out_name}" \
      --output_path_performances="results/performance_2mb_hm450/${out_name}" \
      --output_path_selected_features="results/sel_feats_2mb_hm450/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_2mb_hm450/${out_name}"
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

echo "Running HM450 2Mb-window no-clustering validation for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_missing=$(awk -F'\t' '$2=="missing_prereq"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_missing missing prerequisites, $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/hm450_2mb_validation/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
