"""The CpGs behind the saturation curve: the catalogue itself, not just its size.

Same inputs, same filters and same cohort order as saturation_plot.py, so the cumulative counts it
prints reproduce that run's saturation_counts.tsv (it checks itself against that file when present).
Pass exactly the arguments you passed to saturation_plot.py.

Usage: python scripts/saturation_catalogue.py <WORK_DIR> [--out_suffix=SUFFIX] [--first=LABEL=FILE[:N_SAMPLES]]
                                              [--stable[=TAG]] [--cv=all] [--pairs] [COHORT ...]

Writes into the same <WORK_DIR>/saturation_.../ folder as saturation_plot.py:
  catalogue_cpgs.txt       the catalogue, one CpG per line, sorted (TCGA: cg ids; WGBS-only CpGs: chr_pos)
  catalogue_cpgs.tsv.gz    cpg, n_cohorts, first_step, first_cohort, cohorts  (cohorts comma-separated)
  catalogue_steps.tsv      new / cumulative per step, to compare with saturation_counts.tsv
  catalogue_pairs.tsv.gz   with --pairs: gene, cpg, n_cohorts, first_step, first_cohort
"""
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd

work = sys.argv[1]
SC = os.environ.get("METH_SCALE", "beta")
firsts = [a[len("--first="):] for a in sys.argv[2:] if a.startswith("--first=")]
first = firsts[0] if firsts else None
stable = next(((a.split("=", 1)[1] if "=" in a else "corr0.5") for a in sys.argv[2:] if a.startswith("--stable")), None)
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(
    os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
CV_ALL = any(a == "--cv=all" for a in sys.argv[2:])
PAIRS = any(a == "--pairs" for a in sys.argv[2:])
TH = -np.inf if CV_ALL else 0.1                      # saturation_plot.py draws 0.1 and 0.3; 0.1 is the first column
out = (f"{work}/saturation"
       + ("_" + os.path.basename(first.rsplit("=", 1)[1].split(":")[0]).split(".")[0] if first else "")
       + next((a[len("--out_suffix="):] for a in sys.argv[2:] if a.startswith("--out_suffix=")), "")
       + ("_stable" if stable else "") + ("_cvall" if CV_ALL else ""))
os.makedirs(out, exist_ok=True)


def stable_probes(c):
    """gene -> selected probes lying in a stable co-methylated region (as in saturation_plot.py)."""
    S = f"{work}/{c}/stability_{SC}"
    cpg = pd.read_csv(f"{S}/cpg_freq.tsv", sep="\t", usecols=["TxID", "CpGID", "saved_selected"], dtype={"TxID": str})
    cpg = cpg[cpg.saved_selected == 1]
    rm = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 2, 3], names=["TxID", "CpGID", "region"],
                                dtype={0: str})
                    for f in glob.glob(f"{S}/region_map_{stable}_shard*.tsv")]).drop_duplicates(["TxID", "CpGID"])
    rf = pd.read_csv(f"{S}/region_freq_{stable}.tsv", sep="\t", usecols=["TxID", "region", "stable"], dtype={"TxID": str})
    m = cpg.merge(rm, on=["TxID", "CpGID"]).merge(rf[rf.stable], on=["TxID", "region"])
    return m.groupby("TxID").CpGID.apply(set).to_dict()


def load(c):
    st = stable_probes(c) if stable else None
    rows = []
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        g = os.path.basename(f)[:-8]
        v = dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())
        sf = f"{work}/{c}/results_{SC}/sel_feats/{g}_SVR.tsv"
        if not os.path.exists(sf):
            continue
        probes = [l.strip() for l in open(sf) if l.strip()]
        if st is not None:
            probes = [x for x in probes if x in st.get(g, ())]
            if not probes:
                continue
        rows.append((g, float(v["cv_r2"]), probes))
    return rows


data, nsamp = {}, {}
for c in cohorts:
    data[c] = load(c)
    f = f"{work}/{c}/need_gene_exp_{c}.csv"
    nsamp[c] = len(open(f).readline().strip().split(",")) - 1 if os.path.exists(f) else 0
    print(f"{c}: {len(data[c]):,} gene models", flush=True)
FIRST = []
for fa in firsts:
    lab, f = fa.split("=", 1)
    n = 0
    if ":" in f:
        f, n = f.rsplit(":", 1)
    d = pd.read_csv(f, sep="\t")
    data[lab] = [(g, x.cv_r2.iloc[0], list(x.probe)) for g, x in d.groupby("gene")]
    nsamp[lab] = int(n)
    FIRST.append(lab)
    print(f"{lab} (first): {len(data[lab]):,} gene models with >= 1 selected CpG", flush=True)
cohorts = FIRST + cohorts

# the same order as the bar panel: the external cohorts first, then TCGA by decreasing sample size
order = FIRST + [c for c in sorted(cohorts, key=lambda c: -nsamp[c]) if c not in FIRST]

cpgs = {c: {x for g, cv, ps in data[c] if cv >= TH for x in ps} for c in cohorts}
n_cohorts, first_step, first_cohort, members = {}, {}, {}, {}
seen, rows = set(), []
for step, c in enumerate(order, 1):
    new = cpgs[c] - seen
    for x in new:
        first_step[x], first_cohort[x] = step, c
    seen |= cpgs[c]
    rows.append(dict(step=step, cohort=c, n_samples=nsamp[c], cohort_size=len(cpgs[c]),
                     new=len(new), cumulative=len(seen), new_share_of_cohort=len(new) / max(1, len(cpgs[c]))))
for c in cohorts:
    for x in cpgs[c]:
        n_cohorts[x] = n_cohorts.get(x, 0) + 1
        members.setdefault(x, []).append(c)

ids = sorted(seen)
with open(f"{out}/catalogue_cpgs.txt", "w") as fh:
    fh.write("\n".join(ids) + "\n")
pd.DataFrame({"cpg": ids,
              "n_cohorts": [n_cohorts[x] for x in ids],
              "first_step": [first_step[x] for x in ids],
              "first_cohort": [first_cohort[x] for x in ids],
              "cohorts": [",".join(sorted(members[x], key=order.index)) for x in ids]}
             ).to_csv(f"{out}/catalogue_cpgs.tsv.gz", sep="\t", index=False, compression="gzip")
S = pd.DataFrame(rows)
S.to_csv(f"{out}/catalogue_steps.tsv", sep="\t", index=False, float_format="%.6g")
print(S.to_string(index=False))
print(f"\ncatalogue: {len(ids):,} distinct CpGs, {sum(n_cohorts[x] > 1 for x in ids):,} found in more than one cohort")

if PAIRS:
    pr = {}
    for step, c in enumerate(order, 1):
        for g, cv, ps in data[c]:
            if cv < TH:
                continue
            for x in ps:
                k = (g, x)
                if k not in pr:
                    pr[k] = [0, step, c]
                pr[k][0] += 1
    pd.DataFrame([(g, x, v[0], v[1], v[2]) for (g, x), v in sorted(pr.items())],
                 columns=["gene", "cpg", "n_cohorts", "first_step", "first_cohort"]
                 ).to_csv(f"{out}/catalogue_pairs.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"pairs: {len(pr):,} distinct gene-CpG pairs")

# self-check against the plot's own counts, when that run is present
ref = f"{out}/saturation_counts.tsv"
if os.path.exists(ref):
    R = pd.read_csv(ref, sep="\t")
    R = R[(R.catalogue == "probes") & (R.cv_threshold == (-np.inf if CV_ALL else 0.1))]
    m = R.merge(S, on=["step", "cohort"], suffixes=("_plot", "_here"))
    bad = m[(m.new_plot != m.new_here) | (m.cumulative_plot != m.cumulative_here)]
    print("\nself-check against saturation_counts.tsv:",
          "identical at every step" if bad.empty and len(m) == len(S)
          else f"MISMATCH in {len(bad)} of {len(m)} steps\n" + bad[["step", "cohort", "new_plot", "new_here",
                                                                   "cumulative_plot", "cumulative_here"]].to_string(index=False))
print(f"wrote {out}/catalogue_cpgs.txt and catalogue_cpgs.tsv.gz")
