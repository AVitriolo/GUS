#!/bin/bash
# CDState on mCRPC metastases (PAMES purity), k=3..8 x 10 restarts; resumable.
cd /home/Prostate_GAS
exec /home/ieo7293/miniforge3/envs/cdstate/bin/python analyses/scripts/run_cdstate_multi.py \
  --bulk=results/cdstate/M/bulk_cpm_M.csv --purity=results/cdstate/M/purity_pames_M.csv \
  --outdir=results/cdstate/M/runs --prefix=M --kmin=3 --kmax=8 --niter=10 --ncore=6
