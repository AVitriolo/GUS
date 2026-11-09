#!/usr/bin/env bash
# =============================================================================
# run_methyldriver_gimme_w20.sh
# methyldriver + gimme on the 2 Mb per-CpG validation models
# (results/{sel_feats,performance}_2mb_full_kcapped_cv3r5), bulk data.
#   regions : selected CpGs of the same transcript <= 20 bp apart chained into
#             one region, window = first CpG - 10 .. last CpG + 10 (20 bp for a
#             lone CpG); models with cv_r2 > mean cv_r2; transcripts with exactly
#             the same region collapsed into one
#   bigWigs : read in place from the HPC folder mounted (sshfs, read-only) on
#             resources/roadmap; or, with BW_SOURCES=$DATA/encode_bigwig_s3_urls.tsv,
#             streamed from the ENCODE S3 bucket one at a time (not kept)
#   methyldriver : 100 neighbours from PCA (20 PCs) + kNN, skipping neighbours
#             closer than 1 kb; deltas from the combat H5 (cluster_delta_methylation.R)
#   gimme   : scan of the distinct windows (no padding: all >= 20 bp)
# Nothing existing is overwritten: all outputs carry the _w20 tag.
# Resumable: finished steps are skipped.
#
# Usage: nohup ./run_methyldriver_gimme_w20.sh > logs/methyldriver_w20/run.log 2>&1 &
# =============================================================================
set -euo pipefail
cd /home/Prostate_GAS

MD_ENV=/home/Prostate_GAS/.snakemake/conda/4e4f589e39a867dc7bbcf8f76be7cace_       # workflow/envs/r_methyldriver.yml
R_MERGE_CPGS_ENV=/home/Prostate_GAS/.snakemake/conda/771b4c4756ff2168fbfe4a5a96647c76_  # workflow/envs/r_merge_CpGs.yml
GIMME_ENV=/home/Prostate_GAS/.snakemake/conda/751a681c2274ad9626e43e15de88c327_    # workflow/envs/gimme.yml

TAG=2mb_full_kcapped_cv3r5_w20
DATA=data/to_methyldriver_w20
RES=results/methyldriver_2mb_full_cv3r5              # featmatrix + neighbours from the HPC (build_neighbors_w20.R)
GIM=results/gimme_w20
LOGS=logs/methyldriver_w20
INTEGRATED_H5=resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5
URLS="${BW_SOURCES:-$DATA/bigwig_sources_mount.tsv}"   # accession <TAB> mounted path (or S3 url)
COLUMNS=$DATA/bigwig_columns.txt              # <accession>_<target>, columns of the HPC-built matrix
BW_PARALLEL="${BW_PARALLEL:-3}"               # mounted via the HPC login node: keep it low
BW_TMP=resources/beegfs_scratch/bigwig_w20_tmp
COMPARISONS=(N_vs_T T_vs_D T_vs_D_Bone T_vs_D_Liver T_vs_D_Lymph N_vs_D N_vs_D_Bone N_vs_D_Liver N_vs_D_Lymph)
declare -A C1=([N_vs_T]=N [T_vs_D]=T [T_vs_D_Bone]=T [T_vs_D_Liver]=T [T_vs_D_Lymph]=T [N_vs_D]=N [N_vs_D_Bone]=N [N_vs_D_Liver]=N [N_vs_D_Lymph]=N)
declare -A C2=([N_vs_T]=T [T_vs_D]=D [T_vs_D_Bone]=D [T_vs_D_Liver]=D [T_vs_D_Lymph]=D [N_vs_D]=D [N_vs_D_Bone]=D [N_vs_D_Liver]=D [N_vs_D_Lymph]=D)
declare -A SITE=([T_vs_D_Bone]=Bone [T_vs_D_Liver]=Liver [T_vs_D_Lymph]=Lymph_node [N_vs_D_Bone]=Bone [N_vs_D_Liver]=Liver [N_vs_D_Lymph]=Lymph_node)

mkdir -p $DATA $RES $GIM $LOGS
[ -s $URLS ] && [ -s $COLUMNS ] || { echo "missing $URLS or $COLUMNS" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

# ---- 1. regions ----
if [ ! -s $DATA/md_percpg_$TAG.bed ]; then
    log "1. regions"
    $MD_ENV/bin/Rscript analyses/scripts/sel_feats_to_regions_w20.R \
        --sel_feats_dir=results/sel_feats_2mb_full_kcapped_cv3r5 \
        --perf_dir=results/performance_2mb_full_kcapped_cv3r5 \
        --out_dir=$DATA --tag=$TAG --gap=20 --flank=10 > $LOGS/regions.log 2>&1
    cat $LOGS/regions.log
fi

# ---- 2. deltas (background: independent of the bigWigs) ----
deltas() {
    for c in "${COMPARISONS[@]}"; do
        out=$DATA/${c}_$TAG.tsv
        [ -s $out ] && continue
        echo "$c ${C1[$c]} ${C2[$c]} $out ${SITE[$c]:--}"         # comparison cond_1 cond_2 output site_2 ("-" = none;
                                                                  # no trailing blank: xargs -L would join the next line)
    done | xargs -P 3 -L 1 bash -c '
        '"$R_MERGE_CPGS_ENV"'/bin/Rscript workflow/scripts/cluster_delta_methylation.R \
            --input_h5='"$INTEGRATED_H5"' --input_bed='"$DATA/md_percpg_$TAG.bed"' \
            --cond_1=$1 --cond_2=$2 $([ "$4" != - ] && echo --site_2=$4) --output=$3 > '"$LOGS"'/delta_$0.log 2>&1'
}
# ---- 3. gimme (background) ----
run_gimme() {
    [ -s $GIM/scan_$TAG.tsv ] && return 0
    PATH=$GIMME_ENV/bin:$PATH bash workflow/scripts/gimme_scan.sh \
        $DATA/gimme_coords_$TAG.bed GRCh38.p14 $GIM/scan_$TAG.tsv \
        resources/reference_data/hg38.chrom.sizes 20 "" "" > $LOGS/gimme.log 2>&1
}
log "2. deltas and 3. gimme started in the background"
deltas & pid_delta=$!
run_gimme & pid_gimme=$!

# ---- 4. bigWig averages over the distinct windows, streamed ----
FM=$RES/featmatrix_$TAG.rds
if [ ! -s $FM ]; then
    log "4. bigWig averages ($(wc -l < $URLS) files, $BW_PARALLEL at a time)"
    PATH=$MD_ENV/bin:$PATH bash analyses/scripts/bigwig_stream_average.sh \
        $URLS $DATA/coords_all_$TAG.bed $DATA/bigwig_averages $BW_TMP $BW_PARALLEL > $LOGS/bigwig.log 2>&1
    log "5. feature matrix"
    $MD_ENV/bin/Rscript analyses/scripts/assemble_featmatrix_w20.R \
        --bed=$DATA/coords_all_$TAG.bed --columns=$COLUMNS \
        --part_dir=$DATA/bigwig_averages --output=$FM > $LOGS/featmatrix.log 2>&1
fi

wait $pid_delta || { log "deltas failed, see $LOGS/delta_*.log"; exit 1; }
log "deltas done"

# ---- 6. methyldriver ----
if [ ! -s $RES/full_N_vs_D_Lymph_$TAG.tsv ]; then
    log "6. methyldriver"
    D=$(for c in "${COMPARISONS[@]}"; do printf "%s=%s," $c $DATA/${c}_$TAG.tsv; done)
    $MD_ENV/bin/Rscript analyses/scripts/methyldriver_w20.R \
        --featmatrix=$FM --md_regions=$DATA/md_regions_$TAG.tsv --deltas=${D%,} \
        --conversion_table=resources/conversion_table_gencodeV38.tsv \
        --out_dir=$RES --tag=$TAG --n=100 --n_pcs=20 --min_dist=1000 --q=0.05 \
        --num_cores=16 > $LOGS/methyldriver.log 2>&1
fi

wait $pid_gimme || { log "gimme failed, see $LOGS/gimme.log"; exit 1; }
log "done"
