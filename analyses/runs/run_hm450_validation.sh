#!/usr/bin/env bash
# Manual validation run: how does model performance change when the CpG universe is
# downsampled to Illumina HM450 probes before clustering (merge_CpGs)?
#
# Reuses the existing, UNMODIFIED unpack_RSE_NTM.R / merge_CpGs_split_new.R /
# xgb_input_clustered_split.R / model_selection_split_kbest.py scripts with new
# arguments and new "_hm450"-suffixed output paths, so nothing under data/ or
# results/ that already exists is touched. Not wired into any Snakefile - run by hand.
#
# Usage: ./run_hm450_validation.sh [TxID ...]
#   No args -> loops over every TxID currently present in data/corr_NTM/
#   One or more TxIDs -> restricts the run to just those (useful for a quick spot-check
#   before committing to the full ~27,648-TxID run)
#
# Parallelism: PARALLEL_JOBS TxIDs are processed concurrently via `xargs -P` (GNU
# parallel isn't installed on this machine). Each TxID's own model_selection step
# additionally uses N_JOBS threads internally (RandomSearchCV), so the total thread
# budget is roughly PARALLEL_JOBS * N_JOBS - tune both to your core count, e.g.
# PARALLEL_JOBS=8 ./run_hm450_validation.sh (defaults below assume a shared machine).
#
# Environments: none of the named conda envs on this machine actually match what
# their same-named workflow/envs/*.yml spec requires (e.g. the named "r_merge_CpGs"
# env is missing GenomicRanges/HDF5Array/SummarizedExperiment/ComplexHeatmap/FNN/
# igraph entirely). The envs that actually work are the hashed ones Snakemake itself
# built under .snakemake/conda/ via `--use-conda` - each *_ENV below was identified
# by matching its cached .yaml spec byte-for-byte against the corresponding
# workflow/envs/*.yml, then verified to import every package the relevant script
# needs. If you rebuild any of these (e.g. `snakemake --use-conda --conda-cleanup-envs`,
# or editing the .yml), the hash changes - update the path below to match.

set -uo pipefail

R_XGB_INPUT_ENV=/home/Prostate_GAS/.snakemake/conda/1be302004981d9afc5efa19c4056445a_   # matches workflow/envs/r_xgb_input.yml
R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_  # matches workflow/envs/r_merge_CpGs.yml
PY_ML_ENV=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_         # matches workflow/envs/py_ML.yml

CV=3
N_ITER=27
N_JOBS=2                        # per-TxID sklearn RandomSearchCV threads (kept low since PARALLEL_JOBS adds another dimension)
PARALLEL_JOBS="${PARALLEL_JOBS:-13}"  # concurrent TxIDs; override with e.g. `PARALLEL_JOBS=16 ./run_hm450_validation.sh`

COUNTS_TRAIN=resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv
COUNTS_TEST=resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv
TRAIN_SAMPLES=resources/split/train_samples.txt
TEST_SAMPLES=resources/split/test_samples.txt
INTEGRATED_H5=resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5
PROBE_BED=resources/illumina_probes/HM450_hg38_probes.bed
PREFIX=hg38_v38_pc_50000_25_10

# Every PLOT_EVERY successfully-completed TxIDs, redraw the cv_r2 comparison plot
# (HM450-restricted vs. the original full-CpG results/performance for the same
# TxIDs). PLOT_OUTPUT is a single fixed path - each redraw overwrites it, it is
# never versioned/timestamped.
PLOT_EVERY=500
PERF_HM450_DIR=results/performance_hm450
PERF_ORIGINAL_DIR=results/performance
PLOT_OUTPUT=results/plots/hm450_vs_original_cv_r2.png
PLOT_LOCK_DIR=logs/hm450_validation/_plot.lock

mkdir -p data/rse_NTM_hm450 data/corr_NTM_hm450 resources/CpG_data/top_K_hm450 \
         data/clustering/clusters_maps_hm450 results/plots/clustering_hm450 \
         data/xgb_clustered_hm450 \
         results/performance_hm450 results/sel_feats_hm450 results/hyper_hm450 results/shap_hm450 \
         results/models logs/hm450_validation results/plots

STATUS_FILE="logs/hm450_validation/_status_$$.tsv"
: > "$STATUS_FILE"
rmdir "$PLOT_LOCK_DIR" 2>/dev/null || true   # clear a stale lock left by a killed prior run

# Call each env's own interpreter binary directly rather than `conda run` - this
# env was built by a different conda installation/version (25.3.0) than what's on
# this shell's PATH (26.1.1), and `conda run -p ...` chokes across that version gap
# (a "conda: error: argument --json: ignored explicit argument '1'" parse failure).
# Direct invocation sidesteps it entirely and needs no active/matching conda on PATH.
run_r_xgb_input()   { "$R_XGB_INPUT_ENV/bin/Rscript" "$@"; }
run_r_merge_CpGs()  { "$R_MERGE_CPGS_ENV/bin/Rscript" "$@"; }
run_py_ML()         { "$PY_ML_ENV/bin/python" "$@"; }

process_txid() {
  local TxID="$1"
  local base="${PREFIX}_${TxID}"
  local corr_full="data/corr_NTM/${base}"
  local log="logs/hm450_validation/${base}.log"

  if [ ! -f "$corr_full" ]; then
    printf '%s\tno_corr\n' "$TxID" >> "$STATUS_FILE"
    return
  fi

  # Give each TxID its own TMPDIR under the project rather than relying on the
  # shared system /tmp - a concurrent run hit a transient
  # "cannot open ... HDF5Array_dump_names_global_counter" failure from R's
  # HDF5Array namespace load, most likely system /tmp contention/cleanup under
  # many simultaneous R processes; a private, pre-created TMPDIR avoids it.
  local txid_tmp="logs/hm450_validation/tmp/${base}"
  mkdir -p "$txid_tmp"
  export TMPDIR="$PWD/$txid_tmp"

  (
    echo "=== $TxID ==="

    # 1. derive the probe-restricted CpG list + subset the existing correlation matrix
    run_r_merge_CpGs analyses/scripts/derive_probe_subset_NTM.R \
      --input_path_corr="$corr_full" \
      --probe_bed="$PROBE_BED" \
      --tx_id="$TxID" \
      --output_path_cpg_list="resources/CpG_data/top_K_hm450/${base}.bed" \
      --output_path_corr="data/corr_NTM_hm450/${base}" \
    || exit 1

    # 2. regenerate the per-TxID RSE, restricted to the surviving probe CpGs
    run_r_xgb_input workflow/scripts/unpack_RSE_NTM.R \
      --input_dir="$INTEGRATED_H5" \
      --input_path_filtered_CpGs="resources/CpG_data/top_K_hm450/${base}.bed" \
      --output_path_rse="data/rse_NTM_hm450/${base}" \
    || exit 1

    # 3. re-cluster
    run_r_merge_CpGs workflow/scripts/merge_CpGs_split_new.R \
      --input_path="data/corr_NTM_hm450/${base}" \
      --input_path_rse="data/rse_NTM_hm450/${base}" \
      --input_path_counts="$COUNTS_TRAIN" \
      --train_samples="$TRAIN_SAMPLES" \
      --output_path_plot="results/plots/clustering_hm450/${base}_heatmap.pdf" \
      --output_path_map="data/clustering/clusters_maps_hm450/${base}_clusters_map.txt" \
      --C=2 --j=1 \
    || exit 1

    # 4. rebuild the clustered feature matrix
    run_r_merge_CpGs workflow/scripts/xgb_input_clustered_split.R \
      --input_path_rse="data/rse_NTM_hm450/${base}" \
      --input_path_mapping="data/clustering/clusters_maps_hm450/${base}_clusters_map.txt" \
      --input_path_counts_train="$COUNTS_TRAIN" \
      --input_path_counts_test="$COUNTS_TEST" \
      --train_samples="$TRAIN_SAMPLES" \
      --test_samples="$TEST_SAMPLES" \
      --output_path_train="data/xgb_clustered_hm450/${base}_xgb_input_clustered_train.txt" \
      --output_path_test="data/xgb_clustered_hm450/${base}_xgb_input_clustered_test.txt" \
    || exit 1

    # 5. retrain/evaluate the model
    local out_name="${PREFIX}_${CV}_${TxID}"
    PYTHONPATH=. run_py_ML analyses/scripts/model_selection_split_kbest.py \
      --input_path_dataset_train="data/xgb_clustered_hm450/${base}_xgb_input_clustered_train.txt" \
      --input_path_dataset_test="data/xgb_clustered_hm450/${base}_xgb_input_clustered_test.txt" \
      --input_path_models=config/MODELS \
      --selected_model=SVR \
      --cv="$CV" \
      --optimizer=RandomSearchCV \
      --n_iter="$N_ITER" \
      --n_jobs="$N_JOBS" \
      --output_dir_models=results/models \
      --output_path_hyperparams="results/hyper_hm450/${out_name}" \
      --output_path_performances="results/performance_hm450/${out_name}" \
      --output_path_selected_features="results/sel_feats_hm450/${out_name}" \
      --tx_id="$TxID" \
      --output_path_shap="results/shap_hm450/${out_name}"
  ) > "$log" 2>&1
  local status=$?
  rm -rf "$txid_tmp"

  if [ $status -eq 1 ] && grep -q "skipping (need >= 2" "$log"; then
    printf '%s\tskipped\n' "$TxID" >> "$STATUS_FILE"
  elif [ $status -ne 0 ]; then
    printf '%s\tfailed\n' "$TxID" >> "$STATUS_FILE"
  else
    printf '%s\tok\n' "$TxID" >> "$STATUS_FILE"

    # Every PLOT_EVERY completions, redraw the cv_r2 comparison plot (fixed
    # output path, overwritten). mkdir is atomic, so this lock lets only one
    # concurrent worker actually redraw when several cross the threshold at once.
    local n_ok_now
    n_ok_now=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
    if [ $((n_ok_now % PLOT_EVERY)) -eq 0 ] && mkdir "$PLOT_LOCK_DIR" 2>/dev/null; then
      "$PY_ML_ENV/bin/python" analyses/scripts/plot_hm450_vs_original_cv_r2.py \
        --performance_dir_hm450="$PERF_HM450_DIR" \
        --performance_dir_original="$PERF_ORIGINAL_DIR" \
        --prefix="$PREFIX" \
        --cv="$CV" \
        --output_path="$PLOT_OUTPUT" \
        >> logs/hm450_validation/_plot.log 2>&1
      rmdir "$PLOT_LOCK_DIR"
    fi
  fi
}
export -f process_txid run_r_xgb_input run_r_merge_CpGs run_py_ML
export R_XGB_INPUT_ENV R_MERGE_CPGS_ENV PY_ML_ENV CV N_ITER N_JOBS \
       COUNTS_TRAIN COUNTS_TEST TRAIN_SAMPLES TEST_SAMPLES INTEGRATED_H5 PROBE_BED PREFIX STATUS_FILE \
       PLOT_EVERY PERF_HM450_DIR PERF_ORIGINAL_DIR PLOT_OUTPUT PLOT_LOCK_DIR

if [ "$#" -gt 0 ]; then
  TxIDs=("$@")
else
  TxIDs=()
  for f in data/corr_NTM/${PREFIX}_*; do
    TxIDs+=("$(basename "$f" | sed "s/^${PREFIX}_//")")
  done
fi

echo "Running HM450 validation for ${#TxIDs[@]} TxID(s), ${PARALLEL_JOBS} concurrent"

printf '%s\n' "${TxIDs[@]}" | xargs -P "$PARALLEL_JOBS" -I{} bash -c 'process_txid "$@"' _ {}

n_ok=$(awk -F'\t' '$2=="ok"' "$STATUS_FILE" | wc -l)
n_skipped=$(awk -F'\t' '$2=="skipped" || $2=="no_corr"' "$STATUS_FILE" | wc -l)
n_failed=$(awk -F'\t' '$2=="failed"' "$STATUS_FILE" | wc -l)

echo "Done: $n_ok completed, $n_skipped skipped (too few probe CpGs / no corr file), $n_failed failed"
echo "Per-TxID status: $STATUS_FILE"
if [ "$n_failed" -gt 0 ]; then
  echo "Failed TxIDs (see logs/hm450_validation/<base>.log):"
  awk -F'\t' '$2=="failed" {print "  "$1}' "$STATUS_FILE"
fi
