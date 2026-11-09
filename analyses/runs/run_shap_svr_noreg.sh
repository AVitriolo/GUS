#!/usr/bin/env bash
# SHAP of the original SVR genestd models (refit with the saved hyperparameters), KernelExplainer without l1
# regularisation, for the 9,764 genes with >= 1 stable region; then the SHAP vs stability plots.
# Resumable: genes already complete in a shard file are skipped (keep the shard count at 24).
# Usage: setsid nohup ./run_shap_svr_noreg.sh > logs/shap_svr_noreg/run.log 2>&1 &
set -uo pipefail
cd /home/Prostate_GAS
PY=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_/bin/python     # workflow/envs/py_ML.yml
O=results/stability_selection_2mb_v2_beta_genestd_kcapped; S=$O/shap_svr_noreg; L=logs/shap_svr_noreg; N=24
mkdir -p $L
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for i in $(seq 0 $((N - 1))); do
  nice -n 10 $PY analyses/scripts/shap_svr_refit.py --genes $S/genes.txt \
    --sel_dir results/sel_feats_2mb_v2_beta_genestd_kcapped_cv3r5 --hyper_dir results/hyper_2mb_v2_beta_genestd_kcapped_cv3r5 \
    --stab_dir $O --cand_dir resources/CpG_data/top_K_2Mb_merged_v2 \
    --h5 resources/integrated_wgbs_merged_v2_hdf5/integrated_wgbs_combat_beta.h5 \
    --positions results/test_replication_2mb_v2_beta_genestd_kcapped/h5_cpg_positions.tsv \
    --coldata resources/integrated_wgbs_merged_v2_hdf5/meth_colData.tsv \
    --counts_train resources/RNA/genelevel_standard/genelevel_counts_NTM_train_combat_25_10.tsv \
    --shard $i/$N --out_dir $S >> $L/shard$i.log 2>&1 &
done
wait
grep -l Traceback $L/shard*.log && echo "a shard failed"
$PY analyses/scripts/plot_shap_vs_stability.py --shap_dir $S --stab_dir $O --out_prefix $S/shap_vs_stability
$PY analyses/scripts/plot_shap_distribution.py --shap_dir $S
echo "done"
