#!/usr/bin/env bash
# =============================================================================
# run_methyldriver_stable_genestd.sh
# MethylDriver on the STABLE selected CpGs of the genestd gene models (2Mb v2 beta, k-capped), all 12 comparisons.
#   features : each model's selected CpGs in a stable co-methylated region (r >= 0.5, <= 300 bp, PFER <= 1;
#              stable_sel_feats.py -> results/sel_feats_2mb_v2_beta_genestd_kcapped_stable_corr0.5)
#   regions  : as run_methyldriver_gimme_w20.sh (CpGs <= 20 bp chained, +-10 bp), but ALL models kept
#              (--keep=all: stability is the filter, not cv_r2 > mean)
#   deltas   : cluster_delta_methylation.R on the v2 integration (integrated_wgbs_combat_beta.h5: ComBat BETA, so the
#              deltas are beta differences, not M as in the v1 w20 run); 9 comparisons + metastatic site vs site
#   bigWigs  : read in place from the HPC folder mounted on resources/roadmap (sshfs, read-only), as the w20 run;
#              or BW_SOURCES=<accession TAB url> to stream from ENCODE S3
#   methyldriver : 100 neighbours from PCA (20 PCs) + kNN, skipping neighbours closer than 1 kb (methyldriver_w20.R)
# Nothing existing is overwritten (tag below). Resumable: finished steps are skipped.
# Usage: nohup ./run_methyldriver_stable_genestd.sh > logs/methyldriver_stable_genestd/run.log 2>&1 &
# =============================================================================
set -euo pipefail
cd /home/Prostate_GAS

MD_ENV=/home/Prostate_GAS/.snakemake/conda/4e4f589e39a867dc7bbcf8f76be7cace_       # workflow/envs/r_methyldriver.yml
R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_  # workflow/envs/r_merge_CpGs.yml
PY=/home/Prostate_GAS/.snakemake/conda/018af87306ec3aaba491a8b6ff85ddc0_/bin/python     # workflow/envs/py_ML.yml

TAG=2mb_v2_beta_genestd_stable_corr0.5_w20
STAB=results/stability_selection_2mb_v2_beta_genestd_kcapped
SEL=results/sel_feats_2mb_v2_beta_genestd_kcapped_stable_corr0.5
PERF=results/performance_2mb_v2_beta_genestd_kcapped_cv3r5
DATA=data/to_methyldriver_w20
RES=results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5
LOGS=logs/methyldriver_stable_genestd
INTEGRATED_H5=resources/integrated_wgbs_merged_v2_hdf5/integrated_wgbs_combat_beta.h5
URLS="${BW_SOURCES:-$DATA/bigwig_sources_mount.tsv}"   # accession <TAB> mounted path (or S3 url)
COLUMNS=$DATA/bigwig_columns.txt
BW_PARALLEL="${BW_PARALLEL:-3}"                        # mounted via the HPC login node: keep it low
BW_TMP=resources/beegfs_scratch/bigwig_stable_tmp
# comparison -> cond_1 cond_2 site_1 site_2 ("-" = none)
COMPARISONS=(
  "N_vs_T N T - -"            "T_vs_D T D - -"            "N_vs_D N D - -"
  "T_vs_D_Bone T D - Bone"    "T_vs_D_Liver T D - Liver"  "T_vs_D_Lymph T D - Lymph_node"
  "N_vs_D_Bone N D - Bone"    "N_vs_D_Liver N D - Liver"  "N_vs_D_Lymph N D - Lymph_node"
  "D_Bone_vs_D_Liver D D Bone Liver" "D_Bone_vs_D_Lymph D D Bone Lymph_node" "D_Lymph_vs_D_Liver D D Lymph_node Liver"
)
mkdir -p $DATA $RES $LOGS
log() { echo "[$(date '+%F %T')] $*"; }

# ---- 1. stable features and regions ----
if [ ! -s $DATA/md_percpg_$TAG.bed ]; then
    log "1. stable features and regions"
    $PY analyses/scripts/stable_sel_feats.py --stab_dir $STAB --out_dir $SEL
    $MD_ENV/bin/Rscript analyses/scripts/sel_feats_to_regions_w20.R --sel_feats_dir=$SEL --perf_dir=$PERF \
        --out_dir=$DATA --tag=$TAG --gap=20 --flank=10 --keep=all > $LOGS/regions.log 2>&1
    cat $LOGS/regions.log
fi

# ---- 2. deltas, 12 comparisons (background, 3 at a time) ----
deltas() {
    for c in "${COMPARISONS[@]}"; do
        read -r name c1 c2 s1 s2 <<< "$c"
        [ -s $DATA/${name}_$TAG.tsv ] && continue
        echo "$name $c1 $c2 $s1 $s2"
    done | xargs -r -P 3 -L 1 bash -c '
        '"$R_MERGE_CPGS_ENV"'/bin/Rscript workflow/scripts/cluster_delta_methylation.R \
            --input_h5='"$INTEGRATED_H5"' --input_bed='"$DATA/md_percpg_$TAG.bed"' --cond_1=$1 --cond_2=$2 \
            $([ "$3" != - ] && echo --site_1=$3) $([ "$4" != - ] && echo --site_2=$4) \
            --output='"$DATA"'/$0_'"$TAG"'.tsv > '"$LOGS"'/delta_$0.log 2>&1'
}
log "2. deltas started in the background"
deltas & pid_delta=$!

# ---- 3. bigWig averages over the distinct windows, then the feature matrix ----
FM=$RES/featmatrix_$TAG.rds
if [ ! -s $FM ]; then
    first=$(awk -F'\t' 'NR==1{print $2}' "$URLS")
    if [[ "$first" != http* ]] && [ ! -s "$first" ]; then
        wait $pid_delta || true
        log "STOP: bigWigs not reachable ($first). Mount the HPC bigWig folder on resources/roadmap (sshfs, read-only)"
        log "      or set BW_SOURCES=<accession TAB url> for S3 streaming, then rerun: finished steps are skipped."
        exit 1
    fi
    log "3. bigWig averages ($(wc -l < $URLS) files, $BW_PARALLEL at a time)"
    PATH=$MD_ENV/bin:$PATH bash analyses/scripts/bigwig_stream_average.sh \
        $URLS $DATA/coords_all_$TAG.bed $DATA/bigwig_averages_$TAG $BW_TMP $BW_PARALLEL > $LOGS/bigwig.log 2>&1
    log "   feature matrix"
    $MD_ENV/bin/Rscript analyses/scripts/assemble_featmatrix_w20.R --bed=$DATA/coords_all_$TAG.bed --columns=$COLUMNS \
        --part_dir=$DATA/bigwig_averages_$TAG --output=$FM > $LOGS/featmatrix.log 2>&1
fi

wait $pid_delta || { log "deltas failed, see $LOGS/delta_*.log"; exit 1; }
log "deltas done"

# ---- 4. methyldriver, all comparisons ----
if [ ! -s $RES/full_D_Lymph_vs_D_Liver_$TAG.tsv ]; then
    log "4. methyldriver"
    D=$(for c in "${COMPARISONS[@]}"; do read -r name _ <<< "$c"; printf "%s=%s," $name $DATA/${name}_$TAG.tsv; done)
    $MD_ENV/bin/Rscript analyses/scripts/methyldriver_w20.R \
        --featmatrix=$FM --md_regions=$DATA/md_regions_$TAG.tsv --deltas=${D%,} \
        --conversion_table=resources/conversion_table_gencodeV38.tsv \
        --out_dir=$RES --tag=$TAG --n=100 --n_pcs=20 --min_dist=1000 --q=0.05 \
        --num_cores=16 > $LOGS/methyldriver.log 2>&1
fi
log "done -> $RES"
