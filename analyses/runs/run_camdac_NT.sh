#!/usr/bin/env bash
# Manual run: CAMDAC-style copy-number correction of CPGEA Tumour betas at the
# GISTIC focal peaks ("CPGEA focal lesions" -> Peak Limits). Normal and
# metastasis samples are left untouched.
#
# Not wired into any Snakefile - run by hand. Writes a NEW h5
# (integrated_wgbs_cna.h5, a copy of integrated_wgbs.h5 with only the in-peak
# Tumour rows rewritten, beta and M), so nothing existing is modified.
#
# Usage: ./run_camdac_NT.sh                 (PUREE purity -> integrated_wgbs_cna.h5)
#        PURITY=pames ./run_camdac_NT.sh    (PAMES purity -> integrated_wgbs_cna_pames.h5)
set -euo pipefail

R_ENV=/home/ieo7293/miniforge3/envs/r_pames/bin/Rscript   # has readxl, GenomicRanges, rhdf5

PURITY=${PURITY:-puree}
case $PURITY in
    puree) PURITY_TSV=resources/purity/puree_purity_T.tsv; SFX=""        ;;
    pames) PURITY_TSV=resources/purity/pames_purity_NT.tsv; SFX="_pames" ;;
    *) echo "PURITY must be puree or pames" >&2; exit 1 ;;
esac

H5_DIR=resources/integrated_wgbs_hdf5
PEAKS=resources/CNAs_NT/cna_peaks_NT.tsv
OUT_H5=$H5_DIR/integrated_wgbs_cna${SFX}.h5
SUMMARY=results/cna_correction/camdac_NT${SFX}_summary.tsv
MIN_CN=0.1

mkdir -p logs/cna_correction results/cna_correction

echo "[$(date)] preparing peaks"
$R_ENV analyses/scripts/prepare_cna_peaks_NT.R \
    --input_xlsx=resources/CNAs_NT/41586_2020_2135_MOESM6_ESM.xlsx \
    --sheet="CPGEA focal lesions" \
    --output_path=$PEAKS 2> logs/cna_correction/prepare_cna_peaks_NT.log

echo "[$(date)] copying integrated h5"
cp $H5_DIR/integrated_wgbs.h5 $OUT_H5
cp $H5_DIR/integrated_wgbs_colData.rds   $H5_DIR/integrated_wgbs_cna${SFX}_colData.rds
cp $H5_DIR/integrated_wgbs_rowRanges.rds $H5_DIR/integrated_wgbs_cna${SFX}_rowRanges.rds

echo "[$(date)] correcting"
$R_ENV analyses/scripts/camdac_correct_NT.R \
    --h5=$OUT_H5 \
    --rowranges=$H5_DIR/integrated_wgbs_rowRanges.rds \
    --coldata=$H5_DIR/integrated_wgbs_colData.rds \
    --peaks=$PEAKS \
    --purity=$PURITY_TSV \
    --min_cn=$MIN_CN \
    --output_summary=$SUMMARY 2> logs/cna_correction/camdac_correct_NT${SFX}.log

echo "[$(date)] done -> $OUT_H5, $SUMMARY"
