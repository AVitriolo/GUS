#!/usr/bin/env bash
# Manual run: CAMDAC purification of every CPGEA Tumour, genome-wide
# (n_t from the GISTIC focal peaks, 2 elsewhere; b_n = matched Normal).
# Normal and metastasis columns are copied unchanged.
#
# Not wired into any Snakefile - run by hand. Writes a NEW h5, nothing
# existing is modified. Needs resources/CNAs_NT/cna_peaks_NT.tsv
# (made by run_camdac_NT.sh).
#
# Usage: ./run_purify_NT.sh                 (PAMES purity)
#        PURITY=puree ./run_purify_NT.sh    (PUREE purity)
set -euo pipefail

R_ENV=/home/ieo7293/miniforge3/envs/r_pames/bin/Rscript   # has GenomicRanges, rhdf5

PURITY=${PURITY:-pames}
case $PURITY in
    pames) PURITY_TSV=resources/purity/pames_purity_NT.tsv ;;
    puree) PURITY_TSV=resources/purity/puree_purity_T.tsv  ;;
    *) echo "PURITY must be pames or puree" >&2; exit 1 ;;
esac

H5_DIR=resources/integrated_wgbs_hdf5
OUT=$H5_DIR/integrated_wgbs_purified_${PURITY}
SUMMARY=results/cna_correction/purify_NT_${PURITY}_summary.tsv

mkdir -p logs/cna_correction results/cna_correction

cp $H5_DIR/integrated_wgbs_colData.rds   ${OUT}_colData.rds
cp $H5_DIR/integrated_wgbs_rowRanges.rds ${OUT}_rowRanges.rds

echo "[$(date)] purifying"
$R_ENV analyses/scripts/purify_tumour_NT.R \
    --in_h5=$H5_DIR/integrated_wgbs.h5 \
    --out_h5=${OUT}.h5 \
    --rowranges=$H5_DIR/integrated_wgbs_rowRanges.rds \
    --coldata=$H5_DIR/integrated_wgbs_colData.rds \
    --peaks=resources/CNAs_NT/cna_peaks_NT.tsv \
    --purity=$PURITY_TSV \
    --min_cn=0.1 \
    --output_summary=$SUMMARY 2> logs/cna_correction/purify_NT_${PURITY}.log

echo "[$(date)] done -> ${OUT}.h5, $SUMMARY"
