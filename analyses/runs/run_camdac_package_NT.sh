#!/usr/bin/env bash
# Manual run: CAMDAC package (v0.2.0) purification of every CPGEA Tumour,
# via analyses/scripts/run_camdac_NT.R (meth files from the CPGEA source
# counts, genome-wide CNA segments from the GISTIC peaks, PAMES purity,
# matched Normal as infiltrates proxy).
#
# Not wired into any Snakefile - run by hand. Tumours whose pure.csv.gz
# already exists are skipped, so the script can be re-launched after a crash.
#
# Resources: one pair peaks at ~20 GB RAM and uses ~4-5 threads (data.table),
# so the default is N_JOBS=6 x THREADS=5 = 30 cores, ~120 GB RAM.
#
# Usage: ./run_camdac_package_NT.sh
#        N_JOBS=4 THREADS=5 ./run_camdac_package_NT.sh
set -euo pipefail

R_ENV=/home/ieo7293/miniforge3/envs/r_camdac/bin/Rscript
N_JOBS=${N_JOBS:-6}
THREADS=${THREADS:-5}
OUTDIR=results/camdac_NT
LOGDIR=logs/cna_correction/camdac_package

mkdir -p $OUTDIR $LOGDIR

# every CPGEA Tumour in the integrated colData
TUMOURS=$($R_ENV -e 'cd <- readRDS("resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds");
                     cat(rownames(cd)[cd$condition == "Tumour"], sep = "\n")' 2>/dev/null)

echo "[$(date)] $(echo "$TUMOURS" | wc -l) tumours, $N_JOBS jobs x $THREADS threads"

# round-robin shards, one process per shard
for i in $(seq 0 $((N_JOBS - 1))); do
    SHARD=$(echo "$TUMOURS" | awk -v n=$N_JOBS -v i=$i '(NR - 1) % n == i' | paste -sd,)
    R_DATATABLE_NUM_THREADS=$THREADS OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=1 \
    $R_ENV analyses/scripts/run_camdac_NT.R \
        --source_se=resources/cpgea_wgbs_with_coverage_hg38/se.rds \
        --coldata=resources/integrated_wgbs_hdf5/integrated_wgbs_colData.rds \
        --peaks=resources/CNAs_NT/cna_peaks_NT.tsv \
        --purity=resources/purity/pames_purity_NT.tsv \
        --chrom_sizes=resources/reference_data/hg38.chrom.clean.sizes \
        --outdir=$OUTDIR \
        --n_cores=1 \
        --tumours=$SHARD > $LOGDIR/shard_$i.log 2>&1 &
done
wait

echo "[$(date)] done: $(find $OUTDIR -name '*.pure.csv.gz' | wc -l) pure files in $OUTDIR"
