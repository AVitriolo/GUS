"""Write the TCGA core probes (selected in EVERY given cohort) and the union (selected in >= 1 cohort), per threshold.

A probe is selected in a cohort when a gene model with cv R2 >= threshold selected it ("any" = every model).
Usage: python scripts/core_probes.py <WORK_DIR> [--out_suffix=SUFFIX] [COHORT ...]     (default: every cohort with results)
Writes <WORK_DIR>/core_probes<SUFFIX>/{core,union}_<threshold>.txt and counts.tsv
"""
import sys, os, glob
import pandas as pd

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta")
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
out = f"{work}/core_probes" + next((a[len("--out_suffix="):] for a in sys.argv[2:] if a.startswith("--out_suffix=")), ""); os.makedirs(out, exist_ok=True)

models = {}
for c in cohorts:
    rows = []
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        sf = f"{work}/{c}/results_{SC}/sel_feats/{os.path.basename(f)}"
        if not os.path.exists(sf): continue
        v = dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())
        rows.append((float(v["cv_r2"]), {l.strip() for l in open(sf) if l.strip()}))
    models[c] = rows; print(f"{c}: {len(rows):,} gene models", flush=True)
res = []
for th in (None, 0.1, 0.3):
    lab = "any" if th is None else str(th)
    sel = [set().union(*(s for cv, s in models[c] if th is None or cv >= th)) for c in cohorts]
    core, union = set.intersection(*sel), set().union(*sel)
    for name, s in (("core", core), ("union", union)):
        with open(f"{out}/{name}_{lab}.txt", "w") as fh: fh.write("\n".join(sorted(s)) + "\n")
    res.append(dict(cv_threshold=lab, n_cohorts=len(cohorts), core=len(core), union=len(union)))
R = pd.DataFrame(res); R.to_csv(f"{out}/counts.tsv", sep="\t", index=False); print(R.to_string(index=False))
print(f"cohorts: {' '.join(cohorts)}\nwrote {out}/")
