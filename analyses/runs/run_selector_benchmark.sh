#!/usr/bin/env bash
# Feature-selector benchmark (genestd 2Mb v2 beta, k-capped run): SelectKBest vs SuSiE vs Boruta (elastic net dropped: ~30 min/gene tuning)
# on 306 genes (34 per CV R2 bin), no covariate adjustment.
#   real run: full train + condition-stratified half-samples (100 for SKB, first 20 for SuSiE, none for Boruta)
#   null run: 63 of those genes (7 per bin), train expression shuffled within condition, full train only
# Step 1 benchmark_selectors.py --skip_enet (SKB + shared inputs), step 2 benchmark_susie_boruta.R, per shard.
# Evaluation / summary: benchmark_evaluate.py, benchmark_summary.R.
#
# Usage: ./run_selector_benchmark.sh [n_shards_real=24] [n_shards_null=8]

set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

PY=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_/bin/python   # workflow/envs/py_ML.yml
RS=/home/Prostate_GAS/.snakemake/conda/r_stability/bin/Rscript                         # analyses/envs/r_stability.yml
H5=resources/integrated_wgbs_merged_v2_hdf5/integrated_wgbs_combat_beta.h5
B=results/selector_benchmark_2mb_v2_beta_genestd_kcapped
L=logs/selector_benchmark
N_REAL="${1:-24}"
N_NULL="${2:-8}"
mkdir -p "$B/real" "$B/null" "$L"
cp "$B/bench_train_rows.txt" "$B/bench_test_rows.txt" "$B/real/"
cp "$B/bench_train_rows.txt" "$B/bench_test_rows.txt" "$B/null/"

COMMON=(--sel_dir results/sel_feats_2mb_v2_beta_genestd_kcapped_cv3r5
        --cand_dir resources/CpG_data/top_K_2Mb_merged_v2 --h5 "$H5"
        --positions results/test_replication_2mb_v2_beta_genestd_kcapped/h5_cpg_positions.tsv
        --coldata resources/integrated_wgbs_merged_v2_hdf5/meth_colData.tsv
        --counts_train resources/RNA/genelevel_standard/genelevel_counts_NTM_train_combat_25_10.tsv
        --counts_test resources/RNA/genelevel_standard/genelevel_counts_NTM_test_combat_25_10.tsv)

run_shard() {   # run_shard <real|null> <i> <n>
  local kind=$1 i=$2 n=$3 extra=()
  [ "$kind" = null ] && extra=(--skip_enet --permute --n_pairs 0 --genes "$B/benchmark_genes_null.txt") \
                     || extra=(--skip_enet --n_pairs 50 --genes "$B/benchmark_genes.txt")
  nice -n 10 "$PY" analyses/scripts/benchmark_selectors.py "${COMMON[@]}" "${extra[@]}" \
      --shard "$i/$n" --out_dir "$B/$kind" > "$L/${kind}_py_shard$i.log" 2>&1 &&
  nice -n 10 "$RS" analyses/scripts/benchmark_susie_boruta.R "$B/$kind" "$i" "$H5" 0 1 20 \
      > "$L/${kind}_R_shard$i.log" 2>&1
  echo "$kind shard $i exit $?"
}

# null genes: 7 per CV R2 bin, from the benchmark genes (benchmark_genes.txt is ordered by bin, 34 each)
awk '(NR - 1) % 34 < 7' "$B/benchmark_genes.txt" > "$B/benchmark_genes_null.txt"

for i in $(seq 0 $((N_REAL - 1))); do run_shard real "$i" "$N_REAL" & done
for i in $(seq 0 $((N_NULL - 1))); do run_shard null "$i" "$N_NULL" & done
wait
echo "benchmark done"
