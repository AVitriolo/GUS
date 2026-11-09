#!/bin/bash
# Lambda-blend residual check; detached, checkpointed per chromosome (resumable).
cd /home/Prostate_GAS
exec /home/ieo7293/miniforge3/envs/r_pames/bin/Rscript analyses/scripts/residual_vs_purity_lambda.R \
  --bed_dir=resources/CpG_data/top_K_2Mb_full --outdir=results/purity_checks/lambda --ncore=6
