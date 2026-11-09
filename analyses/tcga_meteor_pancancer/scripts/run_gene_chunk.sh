#!/usr/bin/env bash
# One array task: build the per-gene tables for one chunk of genes, then fit the models PAR_PER_TASK at a time.
# Model = model_selection_split_kbest_repeated_kcapped.py (SVR, RepeatedKFold CV x N_REPEATS, k capped, no SHAP),
# the same script and settings as the prostate WGBS runs. Genes already done are skipped (rerunning a task is safe).
# Usage: run_gene_chunk.sh <COHORT> <chunk file>
set -uo pipefail
COHORT=$1; CHUNK=$2
source "$(dirname "$0")/../config/pipeline.conf"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
C=$WORK_DIR/$COHORT; R=$C/results_$METH_SCALE; DS=$C/datasets_$METH_SCALE/$(basename "$CHUNK" .txt)
mkdir -p $R/{performance,sel_feats,hyper,shap,models,logs} "$DS"

# retry a python step when it fails while importing packages (BeeGFS hiccup), not on real errors
retry_py() {   # retry_py <log> <command...>
  local log=$1; shift
  for try in 1 2 3; do
    "$@" > "$log" 2>&1; rc=$?
    [ $rc -eq 0 ] && return 0
    grep -qE "FileNotFoundError: .*site-packages|ModuleNotFoundError|ImportError" "$log" || return $rc
    sleep $((30 + RANDOM % 90))
  done
  return 1
}
export -f retry_py

# genes of this chunk not done yet
TODO=$DS/todo.txt; : > $TODO
while read -r g; do [ -n "$g" ] && [ ! -f "$R/performance/${g}_SVR.tsv" ] && echo "$g" >> $TODO; done < "$CHUNK"
echo "[$COHORT $(basename $CHUNK)] $(wc -l < $TODO) of $(wc -l < $CHUNK) genes to do"
[ -s $TODO ] || exit 0

retry_py "$DS/build.log" python "$PIPE_DIR/scripts/03_make_gene_datasets.py" --cohort=$COHORT --data_dir=$C --gene_assoc="$PIPE_DIR/resources/HM450.hg19.gene_association.tsv.gz" \
  --genes_file=$TODO --out_dir=$DS --window_bp=$WINDOW_BP --test_size=$TEST_SIZE --scale=$METH_SCALE || { echo "dataset build failed:"; tail -20 "$DS/build.log"; exit 1; }
cat "$DS/build.log"

fit() {
  g=$1
  [ -f "$DS/${g}_train.tsv" ] || { echo -e "$g\tno_dataset" >> "$DS/status.tsv"; return 0; }
  export PYTHONPATH=. PYTHONUNBUFFERED=1
  cd "$PIPE_DIR/scripts" && retry_py "$R/logs/${g}.log" python workflow/scripts/model_selection_split_kbest_repeated_kcapped.py \
    --input_path_dataset_train="$DS/${g}_train.tsv" --input_path_dataset_test="$DS/${g}_test.tsv" \
    --input_path_models="$PIPE_DIR/config/MODELS" --selected_model=SVR --cv=$CV --n_repeats=$N_REPEATS \
    --optimizer=RandomSearchCV --n_iter=$N_ITER --n_jobs=$N_JOBS --output_dir_models="$R/models" \
    --output_path_hyperparams="$R/hyper/${g}_SVR.tsv" --output_path_performances="$R/performance/${g}_SVR.tsv" \
    --output_path_selected_features="$R/sel_feats/${g}_SVR.tsv" --tx_id="$g" --output_path_shap="$R/shap/${g}_SVR.tsv" --skip_shap \
    && echo -e "$g\tok" >> "$DS/status.tsv" || echo -e "$g\tfailed" >> "$DS/status.tsv"
  rm -f "$R/models/${g}_"*.pkl
  [ "$KEEP_GENE_DATASETS" = 1 ] || rm -f "$DS/${g}_train.tsv" "$DS/${g}_test.tsv"
}
export -f fit; export DS R PIPE_DIR CV N_REPEATS N_ITER N_JOBS KEEP_GENE_DATASETS
xargs -P $PAR_PER_TASK -I{} bash -c 'fit "$@"' _ {} < $TODO
echo "[$COHORT $(basename $CHUNK)] done: $(grep -c ok$ $DS/status.tsv 2>/dev/null) ok, $(grep -c failed$ $DS/status.tsv 2>/dev/null) failed, $(grep -c no_dataset$ $DS/status.tsv 2>/dev/null) without dataset"
