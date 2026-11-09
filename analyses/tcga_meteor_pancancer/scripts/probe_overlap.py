"""Overlap of the functional probes between an external cohort (e.g. the prostate WGBS models) and the TCGA cohorts.

A probe is functional in a cohort when a gene model with cv R2 >= threshold selected it (as in saturation_plot.py).
For each threshold:
  * overlap.tsv: the external cohort's probes (450K IDs only) against the TCGA union (probes selected in >= 1 TCGA
    cohort), against the TCGA intersection (selected in every TCGA cohort), and against each TCGA cohort.
  * unique.tsv: for each cohort (external one included), probes selected by that cohort and by no other cohort,
    as an absolute number and as % of the cohort's own functional probes; computed among TCGA cohorts only, and among
    TCGA + the external cohort.

Usage: python scripts/probe_overlap.py <WORK_DIR> --first=LABEL=FILE [COHORT ...]   (default: every cohort with results)
  FILE: TSV with columns gene, cv_r2, probe (as for saturation_plot.py; non-450K CpG IDs are ignored here)
Writes <WORK_DIR>/probe_overlap/{overlap,unique}.tsv (<WORK_DIR>/probe_overlap_<file name>/ with --first)
"""
import sys, os, glob
import pandas as pd

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta")
first = next((a[len("--first="):] for a in sys.argv[2:] if a.startswith("--first=")), None)
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
TH = (0.1, 0.3)
out = f"{work}/probe_overlap" + ("_" + os.path.basename(first.rsplit("=", 1)[1].split(":")[0]).split(".")[0] if first else "") + next((a[len("--out_suffix="):] for a in sys.argv[2:] if a.startswith("--out_suffix=")), ""); os.makedirs(out, exist_ok=True)

def load(c):   # list of (cv_r2, probes) per gene model
    rows = []
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        sf = f"{work}/{c}/results_{SC}/sel_feats/{os.path.basename(f)}"
        if not os.path.exists(sf): continue
        v = dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())
        rows.append((float(v["cv_r2"]), [l.strip() for l in open(sf) if l.strip()]))
    return rows
data = {c: load(c) for c in cohorts}
for c in cohorts: print(f"{c}: {len(data[c]):,} gene models", flush=True)
lab = None
if first:
    lab, f = first.split("=", 1); f = f.split(":")[0]
    d = pd.read_csv(f, sep="\t"); d = d[d.probe.str.startswith("cg")]
    data[lab] = [(x.cv_r2.iloc[0], list(x.probe)) for _, x in d.groupby("gene")]
    print(f"{lab}: {len(data[lab]):,} gene models with >= 1 selected 450K probe", flush=True)

ov, un = [], []
for th in TH:
    S = {c: {p for cv, ps in data[c] if cv >= th for p in ps} for c in data}
    union = set().union(*(S[c] for c in cohorts)); inter = set.intersection(*(S[c] for c in cohorts))
    if lab:
        w = S[lab]
        ov.append(dict(cv_threshold=th, versus="TCGA union (>= 1 cohort)", n_wgbs=len(w), n_other=len(union), n_shared=len(w & union),
                       pct_of_wgbs=100 * len(w & union) / max(1, len(w)), pct_of_other=100 * len(w & union) / max(1, len(union))))
        ov.append(dict(cv_threshold=th, versus="TCGA intersection (all cohorts)", n_wgbs=len(w), n_other=len(inter), n_shared=len(w & inter),
                       pct_of_wgbs=100 * len(w & inter) / max(1, len(w)), pct_of_other=100 * len(w & inter) / max(1, len(inter))))
        for c in cohorts:
            ov.append(dict(cv_threshold=th, versus=c, n_wgbs=len(w), n_other=len(S[c]), n_shared=len(w & S[c]),
                           pct_of_wgbs=100 * len(w & S[c]) / max(1, len(w)), pct_of_other=100 * len(w & S[c]) / max(1, len(S[c]))))
    for among, group in [("TCGA only", cohorts)] + ([("TCGA + " + lab, cohorts + [lab])] if lab else []):
        for c in group:
            rest = set().union(*(S[o] for o in group if o != c))
            u = len(S[c] - rest)
            un.append(dict(cv_threshold=th, among=among, cohort=c, n_functional=len(S[c]), n_unique=u, pct_unique=100 * u / max(1, len(S[c]))))
    print(f"cv >= {th}: TCGA union {len(union):,}, intersection {len(inter):,}", flush=True)

if ov:
    O = pd.DataFrame(ov); O.to_csv(f"{out}/overlap.tsv", sep="\t", index=False, float_format="%.2f"); print(O.to_string(index=False, float_format="%.2f"))
U = pd.DataFrame(un); U.to_csv(f"{out}/unique.tsv", sep="\t", index=False, float_format="%.2f"); print(U.to_string(index=False, float_format="%.2f"))
print(f"wrote {out}/overlap.tsv, {out}/unique.tsv")
