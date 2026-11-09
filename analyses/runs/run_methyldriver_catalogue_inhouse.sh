#!/usr/bin/env bash
# =============================================================================
# run_methyldriver_catalogue_inhouse.sh
# MethylDriver, in-house gastric cohort metastasis vs tumour (11 matched pairs), on the regions of the stable
# functional-CpG catalogue (saturation_wgbs_genestd_stable_footprint_gbm_stable_cvall: Prostate WGBS + GBM WGBS +
# 15 TCGA cohorts, all stable CpGs, no cv R2 filter; 607,451 entries -> 516,569 regions, catalogue_to_md_regions.R).
#   deltas   : inhouse_delta_methylation.R (>= 5x in >= 80% of pairs; paired, paired t, unpaired), variants
#              all / no_pair8 / no_pair9 (already computed in data/to_methyldriver_w20/inhouse_MvsT_*_<TAG>_*)
#   features : $FM from the HPC (build_neighbors_w20.R, hpc_featmatrix_catalogue), used as is if present; else
#              FULL=1 (default): all 516,569 windows averaged on the HPC (hpc_bigwig_catalogue_all/run_on_hpc.sh),
#              averages copied to $AVG, matrix built from them alone (coords_all_<TAG>.bed);
#              FULL=0: rows of the existing w20 feature matrices where the window is identical, the other 135,527
#              windows averaged on the HPC (hpc_bigwig_catalogue/run_on_hpc.sh)
#   genes    : nearest protein-coding TSS of each region (the catalogue has no gene models)
# Usage: nohup ./run_methyldriver_catalogue_inhouse.sh > logs/methyldriver_catalogue_inhouse/run.log 2>&1 &
# =============================================================================
set -euo pipefail
cd /home/Prostate_GAS
MD_ENV=/home/Prostate_GAS/.snakemake/conda/4e4f589e39a867dc7bbcf8f76be7cace_   # workflow/envs/r_methyldriver.yml
TAG=catalogue_gbm_stable_cvall_w20
DATA=data/to_methyldriver_w20
FULL=${FULL:-1}
AVG=${AVG:-/home/ieo7293/hpc_bigwig_catalogue$([ "$FULL" = 1 ] && echo _all)/averages}
RES=results/methyldriver_${TAG}_inhouse_MvsT
LOGS=logs/methyldriver_catalogue_inhouse
FM=$RES/featmatrix_$TAG.rds
mkdir -p $RES $LOGS
log() { echo "[$(date '+%F %T')] $*"; }

# a matrix built on the HPC with build_neighbors_w20.R (hpc_featmatrix_catalogue), copied to $FM, is used as is
if [ ! -s $FM ]; then
  n=$(ls $AVG/*.tsv.gz 2>/dev/null | wc -l)
  [ "$n" -eq "$(wc -l < $DATA/bigwig_columns.txt)" ] || { log "STOP: no $FM and $n of $(wc -l < $DATA/bigwig_columns.txt) bigWig averages in $AVG"; exit 1; }
fi

if [ ! -s $FM ] && [ "$FULL" = 1 ]; then
  log "feature matrix: all windows from $AVG"
  $MD_ENV/bin/Rscript analyses/scripts/assemble_featmatrix_w20.R --bed=$DATA/coords_all_$TAG.bed --columns=$DATA/bigwig_columns.txt \
      --part_dir=$AVG --output=$FM > $LOGS/featmatrix.log 2>&1
elif [ ! -s $FM ]; then
  log "feature matrix: new windows + existing rows"
  $MD_ENV/bin/Rscript analyses/scripts/assemble_featmatrix_w20.R --bed=$DATA/coords_new_$TAG.bed --columns=$DATA/bigwig_columns.txt \
      --part_dir=$AVG --output=$RES/featmatrix_new_$TAG.rds > $LOGS/featmatrix.log 2>&1
  $MD_ENV/bin/Rscript -e '
    suppressPackageStartupMessages(library(data.table)); D <- "'$DATA'"; TAG <- "'$TAG'"
    need <- unique(unlist(lapply(c("all","no_pair8","no_pair9"), function(v) fread(sprintf("%s/inhouse_MvsT_%s_%s_regions.tsv", D, v, TAG), select = "coordID")$coordID)))
    parts <- list(); left <- need
    for (f in c("'$RES'/featmatrix_new_'$TAG'.rds", "results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5/featmatrix_2mb_v2_beta_genestd_stable_corr0.5_w20.rds",
                "results/methyldriver_2mb_full_cv3r5/featmatrix_2mb_full_kcapped_cv3r5_w20.rds")) {
      x <- readRDS(f); k <- intersect(left, rownames(x)); parts[[f]] <- x[k, , drop = FALSE]; left <- setdiff(left, k); rm(x); invisible(gc()) }
    if (length(left)) stop(length(left), " regions without features")
    cn <- colnames(parts[[1]]); stopifnot(all(vapply(parts, function(p) identical(colnames(p), cn), NA)))
    fm <- do.call(rbind, parts); cat(nrow(fm), "regions x", ncol(fm), "features\n"); saveRDS(fm, "'$FM'")' >> $LOGS/featmatrix.log 2>&1
fi

for v in all no_pair8 no_pair9; do
  P=$DATA/inhouse_MvsT_${v}_$TAG
  [ -s $RES/$v/full_MvsT_unpaired_${v}_$TAG.tsv ] && continue
  $MD_ENV/bin/Rscript analyses/scripts/methyldriver_w20.R --featmatrix=$FM --md_regions=${P}_regions.tsv \
      --deltas=MvsT_paired=${P}_paired.tsv,MvsT_paired_t=${P}_paired_t.tsv,MvsT_unpaired=${P}_unpaired.tsv \
      --conversion_table=resources/conversion_table_gencodeV38.tsv --out_dir=$RES/$v --tag=${v}_$TAG \
      --n=100 --n_pcs=20 --min_dist=1000 --q=0.05 --num_cores=16 > $LOGS/methyldriver_$v.log 2>&1 &
done
wait
grep -h "significant" $LOGS/methyldriver_*.log
log "done -> $RES"
