#!/usr/bin/env bash
# =============================================================================
# run_methyldriver_purebeta_camdac_nocombat.sh
# MethylDriver N vs T and bone vs lymph node on NON-ComBat betas, three versions on identical CpGs/samples/missing cells:
# bulk (v2 raw), CAMDAC, PureBeta (build_purebeta_md_store.R). Regions, feature matrix and neighbours as the earlier
# CAMDAC rerun (results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5_purified_NT), so only the values differ.
# Usage: nohup ./run_methyldriver_purebeta_camdac_nocombat.sh > logs/methyldriver_nocombat/run.log 2>&1 &
# =============================================================================
set -euo pipefail
cd /home/Prostate_GAS
PB_ENV=/home/ieo7293/miniforge3/envs/r_purify
MD_ENV=/home/Prostate_GAS/.snakemake/conda/4e4f589e39a867dc7bbcf8f76be7cace_
R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_
TAG=2mb_v2_beta_genestd_stable_corr0.5_w20
OLD=data/to_methyldriver_w20_purified; BED=$OLD/md_percpg_purified_$TAG.bed
FM=results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5/featmatrix_$TAG.rds
NB=results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5_purified_NT/purified/neighbors_$TAG.rds
D=data/to_methyldriver_w20_nocombat; RES=results/methyldriver_${TAG%_w20}_nocombat; LOGS=logs/methyldriver_nocombat
mkdir -p $D $RES $LOGS
log() { echo "[$(date '+%F %T')] $*"; }

if [ ! -s $D/purebeta.h5 ]; then
    log "1. stores (bulk, CAMDAC, PureBeta)"
    $PB_ENV/bin/Rscript analyses/scripts/build_purebeta_md_store.R --bed=$BED --out_dir=$D --cores=80 2>&1 \
        | grep --line-buffered -v '^  |\|^=\|^$\|Using\|Correcting\|Generating\|PROCESS' > $LOGS/stores.log
    cat $LOGS/stores.log
fi

log "2. deltas"
for v in bulk camdac purebeta; do
    for c in "N_vs_T N T - -" "D_Bone_vs_D_Lymph D D Bone Lymph_node"; do
        read -r name c1 c2 s1 s2 <<< "$c"; o=$D/${name}_${v}_$TAG.tsv; [ -s $o ] && continue
        $R_MERGE_CPGS_ENV/bin/Rscript workflow/scripts/cluster_delta_methylation.R --input_h5=$D/$v.h5 --input_bed=$BED \
            --cond_1=$c1 --cond_2=$c2 $([ $s1 != - ] && echo --site_1=$s1) $([ $s2 != - ] && echo --site_2=$s2) \
            --output=$o > $LOGS/delta_${name}_$v.log 2>&1 &
    done
done
wait; log "   deltas done"

log "3. methyldriver"
for v in bulk camdac purebeta; do
    o=$RES/$v; mkdir -p $o; [ -s $o/full_D_Bone_vs_D_Lymph_$TAG.tsv ] && continue
    cp -n $NB $o/
    $MD_ENV/bin/Rscript analyses/scripts/methyldriver_w20.R --featmatrix=$FM --md_regions=$OLD/md_regions_$TAG.tsv \
        --deltas=N_vs_T=$D/N_vs_T_${v}_$TAG.tsv,D_Bone_vs_D_Lymph=$D/D_Bone_vs_D_Lymph_${v}_$TAG.tsv \
        --conversion_table=resources/conversion_table_gencodeV38.tsv --out_dir=$o --tag=$TAG \
        --n=100 --n_pcs=20 --min_dist=1000 --q=0.05 --num_cores=16 > $o/methyldriver.log 2>&1 &
done
wait; log "done -> $RES"
