#!/usr/bin/env bash
# Does the k-cap generalize better OUT OF COHORT? Trains the WGBS models for the
# HM450-restricted 2Mb candidate set twice per TxID - uncapped
# (model_selection_split_kbest_repeated.py) and k-capped (*_kcapped.py), both SVR,
# RepeatedKFold cv=3/n_repeats=5, RandomSearchCV n_iter=27, same seed - and KEEPS the
# fitted pipelines, so analyses/scripts/apply_wgbs_models_tcga.py can score the exact
# WGBS-trained models on TCGA PRAD (HM450). HM450-restricted because the full-CpG
# models' selected features are ~98% non-array CpGs, i.e. not measurable in TCGA.
#
# Candidates (resources/CpG_data/top_K_2Mb_tcga) = resources/CpG_data/top_K_2Mb further
# restricted to probes present in /home/ieo7293/TCGA_PRAD/PRAD_data_meth.csv: 89,512
# HM450 probes are all-NA in TCGA PRAD (Xena-masked) and the preprocessing dropped
# probes with >20% NA, so without this a model could select CpGs TCGA never measured
# (and the uncapped model, selecting more, would be penalized more for it).
#
# TxID list (resources/TxIDs/tcga_transfer_top50_independent.tsv): greedily, the 50 TxIDs with the most
# TCGA-measured candidates (193-574) sharing NO candidate CpG with each other, one per gene,
# so they are independent tests and the cap (per-fold train size
# 239*2/3 = 159) binds for all of them.
#
# Unlike the validation scripts, the RSE is deleted but the (small) feature matrices
# and the .pkl models are kept - the apply step needs both.
#
# Usage: ./run_tcga_transfer_hm450_2mb_kcap_vs_uncapped.sh resources/TxIDs/tcga_transfer_top50_independent.tsv

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
PARALLEL_JOBS="${PARALLEL_JOBS:-10}"

COUNTS_TRAIN=resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv
COUNTS_TEST=resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv
TRAIN_SAMPLES=resources/split/train_samples.txt
TEST_SAMPLES=resources/split/test_samples.txt
INTEGRATED_H5=resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5
PREFIX=hg38_v38_pc_50000_25_10

OUT=results/tcga_transfer_hm450_2mb
XGB=data/xgb_tcga_transfer_hm450_2mb
RSE=data/rse_tcga_transfer_hm450_2mb
LOGS=logs/tcga_transfer_hm450_2mb

if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  echo "Usage: $0 txid_list.tsv  (first column TxID, header line skipped)" >&2
  exit 1
fi

for cond in uncapped kcapped; do
  mkdir -p "$OUT/$cond"/{performance,hyper,sel_feats,shap,models}
done
mkdir -p "$XGB" "$RSE" "$LOGS" resources/CpG_data/top_K_2Mb_tcga_identity

STATUS_FILE="$LOGS/_status_$$.tsv"
: > "$STATUS_FILE"

run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local cpg_list="resources/CpG_data/top_K_2Mb_tcga/${base}.bed"
  local train_xgb="$XGB/${base}_xgb_input_clustered_train.txt"
  local test_xgb="$XGB/${base}_xgb_input_clustered_test.txt"
  local log="$LOGS/${base}.log"

  if [ ! -f "$cpg_list" ]; then
    printf '%s\tall\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  local txid_tmp="$LOGS/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ==="
    if [ ! -f "$train_xgb" ] || [ ! -f "$test_xgb" ]; then
      run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
        --input_dir="$INTEGRATED_H5" \
        --input_path_filtered_CpGs="$cpg_list" \
        --output_path_rse="$RSE/${base}" \
      || exit 1

      local identity_map="resources/CpG_data/top_K_2Mb_tcga_identity/${base}_identity_map.txt"
      awk -F'\t' '{print $2"\t"$2}' "$cpg_list" > "$identity_map"

      run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
        --input_path_rse="$RSE/${base}" \
        --input_path_mapping="$identity_map" \
        --input_path_counts_train="$COUNTS_TRAIN" \
        --input_path_counts_test="$COUNTS_TEST" \
        --train_samples="$TRAIN_SAMPLES" \
        --test_samples="$TEST_SAMPLES" \
        --output_path_train="$train_xgb" \
        --output_path_test="$test_xgb" \
      || exit 1
      rm -rf "$RSE/${base}"
    else
      echo "[$TxID] reusing existing feature matrix"
    fi
  ) > "$log" 2>&1
  if [ $? -ne 0 ]; then
    printf '%s\tall\tfailed_inputs\n' "$TxID" >> "$STATUS_FILE"
    rm -rf "$txid_tmp"
    return
  fi

  local cond script
  for cond in uncapped kcapped; do
    local out_name="${PREFIX}_${CV}_${TxID}"
    if [ -f "$OUT/$cond/performance/${out_name}" ] && ls "$OUT/$cond/models/${TxID}_"*.pkl >/dev/null 2>&1; then
      printf '%s\t%s\talready_done\n' "$TxID" "$cond" >> "$STATUS_FILE"
      continue
    fi
    if [ "$cond" = kcapped ]; then
      script=workflow/scripts/model_selection_split_kbest_repeated_kcapped.py
    else
      script=analyses/scripts/model_selection_split_kbest_repeated.py
    fi

    PYTHONPATH=. run_py_ML "$script" \
      --input_path_dataset_train="$train_xgb" \
      --input_path_dataset_test="$test_xgb" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --n_repeats="$N_REPEATS" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models="$OUT/$cond/models" \
      --output_path_hyperparams="$OUT/$cond/hyper/${out_name}" \
      --output_path_performances="$OUT/$cond/performance/${out_name}" \
      --output_path_selected_features="$OUT/$cond/sel_feats/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="$OUT/$cond/shap/${out_name}" \
      --skip_shap \
      >> "$log" 2>&1
    if [ $? -ne 0 ]; then
      printf '%s\t%s\tfailed\n' "$TxID" "$cond" >> "$STATUS_FILE"
    else
      printf '%s\t%s\tok\n' "$TxID" "$cond" >> "$STATUS_FILE"
    fi
    rm -f "$OUT/$cond/shap/${out_name}"   # all-NaN placeholder under --skip_shap
  done
  rm -rf "$txid_tmp"
}
export -f process_txid run_r_xgb_input run_r_merge_CpGs run_py_ML
export R_XGB_INPUT_ENV R_MERGE_CPGS_ENV PY_ML_ENV CV N_REPEATS N_ITER N_JOBS \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES INTEGRATED_H5 PREFIX \
       OUT XGB RSE LOGS STATUS_FILE

TxIDs=$(awk -F'\t' 'NR>1 && $1!="" {print $1}' "$1")
echo "WGBS HM450-2Mb uncapped vs k-capped training for $(echo "$TxIDs" | wc -l) TxID(s), ${PARALLEL_JOBS} concurrent - started $(date)"

echo "$TxIDs" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

for s in ok already_done missing_prereq failed_inputs failed; do
  echo "$s: $(awk -F'\t' -v s="$s" '$3==s' "$STATUS_FILE" | wc -l)"
done
echo "Per-TxID status: $STATUS_FILE - finished $(date)"
