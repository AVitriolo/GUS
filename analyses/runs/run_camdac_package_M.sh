#!/usr/bin/env bash
# Manual run: CAMDAC package (v0.2.0) purification of the WCDT metastases
# (liver excluded) via analyses/scripts/run_camdac_M.R: PAMES purity, own or
# site-pooled benign normal, copycat CN converted with purity and ploidy.
# Not wired into any Snakefile - run by hand. Finished metastases are skipped.
#
# Resources: one pair peaks at ~20 GB RAM; N_JOBS=7 -> ~150 GB.
set -euo pipefail
R_ENV=/home/ieo7293/miniforge3/envs/r_camdac/bin/Rscript
N_JOBS=${N_JOBS:-7}
THREADS=${THREADS:-4}
OUTDIR=results/camdac_M
LOGDIR=logs/cna_correction/camdac_M
mkdir -p $OUTDIR $LOGDIR

METS=$($R_ENV -e 'cd <- readRDS("resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds");
                  ex <- readLines("results/camdac_M/excluded_liver.txt");
                  cat(setdiff(rownames(cd)[cd$condition == "metastasis"], ex), sep = "\n")' 2>/dev/null)
echo "[$(date)] $(echo "$METS" | wc -l) metastases, $N_JOBS jobs x $THREADS threads"

# build the shared pooled-normal files once, before the parallel jobs start
$R_ENV analyses/scripts/run_camdac_M.R --mcrpc=resources/mcrpc_wgbs_hg38 \
    --coldata=resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds \
    --exclude=results/camdac_M/excluded_liver.txt --normal_dir=resources/normal_M_wgbs/merged \
    --normal_sites=resources/normal_M_wgbs/normal_sites.tsv --purity=resources/purity/pames_purity_M.tsv \
    --copycat_dir=resources/CNAs_M --copycat_samples=results/camdac_M/copycat_samples.tsv \
    --chrom_sizes=resources/reference_data/hg38.chrom.clean.sizes --outdir=$OUTDIR \
    --tumours=DTB_005,DTB_011,DTB_055 > $LOGDIR/warmup.log 2>&1

for i in $(seq 0 $((N_JOBS - 1))); do
    SHARD=$(echo "$METS" | awk -v n=$N_JOBS -v i=$i '(NR - 1) % n == i' | paste -sd,)
    R_DATATABLE_NUM_THREADS=$THREADS OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=1 \
    $R_ENV analyses/scripts/run_camdac_M.R --mcrpc=resources/mcrpc_wgbs_hg38 \
        --coldata=resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds \
        --exclude=results/camdac_M/excluded_liver.txt --normal_dir=resources/normal_M_wgbs/merged \
        --normal_sites=resources/normal_M_wgbs/normal_sites.tsv --purity=resources/purity/pames_purity_M.tsv \
        --copycat_dir=resources/CNAs_M --copycat_samples=results/camdac_M/copycat_samples.tsv \
        --chrom_sizes=resources/reference_data/hg38.chrom.clean.sizes --outdir=$OUTDIR \
        --tumours=$SHARD > $LOGDIR/shard_$i.log 2>&1 &
done
wait
echo "[$(date)] done: $(find $OUTDIR -name '*.pure.csv.gz' | wc -l) pure files"
