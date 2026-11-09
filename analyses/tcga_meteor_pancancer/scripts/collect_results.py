"""Summarise the results of one or more cohorts: per-gene table + medians/means and a progress line.
Usage: python scripts/collect_results.py <WORK_DIR> BRCA [LUAD ...]   (env METH_SCALE=beta|M, default beta)   -> <WORK_DIR>/<COHORT>/summary_<COHORT>.tsv
"""
import sys, glob, os
import pandas as pd
SC = os.environ.get("METH_SCALE", "beta")
work, cohorts = sys.argv[1], sys.argv[2:]
for c in cohorts:
    rows = []
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        v = dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())
        rows.append(dict(gene=os.path.basename(f)[:-8], **{k: float(v[k]) for k in ("train_r2", "cv_r2", "test_r2") if k in v}))
    n_all = sum(1 for _ in open(f"{work}/{c}/gene_list_{c}.txt")) if os.path.exists(f"{work}/{c}/gene_list_{c}.txt") else float("nan")
    if not rows: print(f"{c}: no results yet"); continue
    d = pd.DataFrame(rows).sort_values("gene"); d.to_csv(f"{work}/{c}/summary_{c}_{SC}.tsv", sep="\t", index=False)
    cl = d[["cv_r2", "test_r2"]].clip(lower=-1)
    print(f"{c}: {len(d):,} / {n_all:,} genes | cv R2 median {d.cv_r2.median():.3f} mean {cl.cv_r2.mean():.3f} | "
          f"test R2 median {d.test_r2.median():.3f} mean {cl.test_r2.mean():.3f} | cv >= 0.1: {100*(d.cv_r2>=0.1).mean():.1f}%")
