#!/usr/bin/env bash
# =============================================================================
# run_methyldriver_stable_inhouse.sh
# MethylDriver, in-house cohort metastasis vs tumour (11 matched pairs), on the
# stable genestd regions of run_methyldriver_stable_genestd.sh.
#   deltas  : analyses/scripts/inhouse_delta_methylation.R (beta, cov >= 5x in
#             >= 80% of pairs, both samples; paired, paired t, unpaired;
#             delta = metastasis - tumour), variants all / no_pair8 / no_pair9
#   regions : only the regions with a delta; neighbours rebuilt among them from
#             the existing chromatin feature matrix (no bigWigs needed)
# Usage: nohup ./run_methyldriver_stable_inhouse.sh > logs/methyldriver_stable_inhouse/run.log 2>&1 &
# =============================================================================
set -euo pipefail
cd /home/Prostate_GAS
MD_ENV=/home/Prostate_GAS/.snakemake/conda/4e4f589e39a867dc7bbcf8f76be7cace_   # workflow/envs/r_methyldriver.yml
T=2mb_v2_beta_genestd_stable_corr0.5_w20
DATA=data/to_methyldriver_w20
FM=results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5/featmatrix_$T.rds
RES=results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5_inhouse_MvsT
LOGS=logs/methyldriver_stable_inhouse
mkdir -p $RES $LOGS
for v in all no_pair8 no_pair9; do
  P=$DATA/inhouse_MvsT_${v}_$T
  [ -s $RES/$v/full_MvsT_unpaired_${v}_$T.tsv ] && continue
  $MD_ENV/bin/Rscript analyses/scripts/methyldriver_w20.R --featmatrix=$FM --md_regions=${P}_regions.tsv \
      --deltas=MvsT_paired=${P}_paired.tsv,MvsT_paired_t=${P}_paired_t.tsv,MvsT_unpaired=${P}_unpaired.tsv \
      --conversion_table=resources/conversion_table_gencodeV38.tsv --out_dir=$RES/$v --tag=${v}_$T \
      --n=100 --n_pcs=20 --min_dist=1000 --q=0.05 --num_cores=16 > $LOGS/methyldriver_$v.log 2>&1 &
done
wait
echo "done -> $RES"
