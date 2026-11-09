#!/bin/bash
# CDState diagnostics batch: per-site, PURPLE-input and CPGEA positive control. Resumable.
cd /home/Prostate_GAS
for d in T_pames M_pames_all M_pames_Bone M_pames_Lymph_node M_purple_all M_purple_Bone M_purple_Lymph_node; do
  echo "=== $d $(date +%T) ==="
  /home/ieo7293/miniforge3/envs/cdstate/bin/python analyses/scripts/run_cdstate_multi.py \
    --bulk=results/cdstate/$d/bulk_cpm.csv --purity=results/cdstate/$d/purity.csv \
    --outdir=results/cdstate/$d/runs --prefix=$d --kmin=3 --kmax=8 --niter=10 --ncore=6 2>&1 | grep -vE "\[INFO\] Error"
done
echo "BATCH FINISHED $(date +%T)"
