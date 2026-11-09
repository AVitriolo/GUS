#!/usr/bin/env bash
# CV-stability ablation: reuses the ALREADY-CACHED 2Mb no-clustering feature
# matrices (data/xgb_2mb_hm450/, the condition at parity with the original
# pipeline) and reruns only the modeling step with
# model_selection_split_kbest_repeated.py (RepeatedKFold instead of a single
# plain KFold - see that script's header and
# workflow/scripts/helpers/OPTIMIZE_HYPERPARAMS_repeated.py for why). No RSE,
# correlation, or clustering work is redone - this is purely a modeling-stage
# comparison against the existing results/performance_2mb_hm450/ (cv=3,
# n_repeats=1, i.e. today's plain KFold) baseline.
#
# Usage: ./run_cv_stability_validation.sh CV N_REPEATS TxID [TxID ...]
#   e.g. ./run_cv_stability_validation.sh 3 5 $(cat /tmp/200_txids.txt)
#        ./run_cv_stability_validation.sh 5 5 $(cat /tmp/200_txids.txt)

set -uo pipefail

PY_ML_ENV=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_  # matches workflow/envs/py_ML.yml

N_ITER=27
N_JOBS=2
PARALLEL_JOBS="${PARALLEL_JOBS:-13}"

PREFIX=hg38_v38_pc_50000_25_10
XGB_DIR=data/xgb_2mb_hm450

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 CV N_REPEATS TxID [TxID ...]" >&2
  exit 1
fi
CV="$1"; shift
N_REPEATS="$1"; shift
TxIDs=("$@")
TAG="cv${CV}r${N_REPEATS}"

PERF_DIR="results/performance_${TAG}"
SELF_DIR="results/sel_feats_${TAG}"
HYPER_DIR="results/hyper_${TAG}"
SHAP_DIR="results/shap_${TAG}"
LOG_DIR="logs/cv_stability_${TAG}"

mkdir -p "$PERF_DIR" "$SELF_DIR" "$HYPER_DIR" "$SHAP_DIR" results/models "$LOG_DIR"

STATUS_FILE="${LOG_DIR}/_status_$$.tsv"
: > "$STATUS_FILE"

run_py_ML() { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local train="${XGB_DIR}/${base}_xgb_input_clustered_train.txt"
  local test="${XGB_DIR}/${base}_xgb_input_clustered_test.txt"
  local log="${LOG_DIR}/${base}.log"

  if [ ! -f "$train" ] || [ ! -f "$test" ]; then
    printf '%s\tmissing_prereq\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  (
    echo "=== $TxID (cv=$CV, n_repeats=$N_REPEATS) ==="
    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest_repeated.py \
      --input_path_dataset_train="$train" \
      --input_path_dataset_test="$test" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --n_repeats="$N_REPEATS" \
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

  if [ $status -ne 0 ]; then
    printf '%s\tfailed\n' "$TxID" >> "$STATUS_FILE"
  else
    printf '%s\tok\n' "$TxID" >> "$STATUS_FILE"
  fi
}
export -f process_txid run_py_ML
export PY_ML_ENV CV N_REPEATS N_ITER N_JOBS PREFIX XGB_DIR STATUS_FILE \
       PERF_DIR SELF_DIR HYPER_DIR SHAP_DIR LOG_DIR

echo "Running CV-stability validation (cv=$CV, n_repeats=$N_REPEATS) for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

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
