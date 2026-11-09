"""Saturation of the catalogue of functional CpGs across TCGA cohorts.

A probe is "functional" in a cohort when a gene model with cv R2 >= threshold selected it. Two catalogues:
  * probes          : distinct 450K probes
  * gene-probe pairs: (gene, probe), i.e. a probe counted once per gene it helps predict
The cumulative size of each catalogue is drawn as cohorts are added, averaged over random cohort orders (band = 5-95%),
and separately in order of decreasing sample size (bars: already in the catalogue vs new).

Usage: python scripts/saturation_plot.py <WORK_DIR> [--out_suffix=SUFFIX] [--first=LABEL=FILE[:N_SAMPLES]] [--stable[=TAG]] [COHORT ...]   (default: every cohort with results)
  --cv=all: no cv R2 threshold (every model counts); one column per catalogue; output folder suffix _cvall.
  --stable: a model contributes only its selected probes lying in a STABLE co-methylated region (stability selection,
           <cohort>/stability_<scale>/: cpg_freq.tsv, region_map_<TAG>_shard*.tsv, region_freq_<TAG>.tsv; TAG default corr0.5);
           models without one are left out. Output folder gets the suffix _stable. Give --first a stable WGBS file too.
  --first: an external cohort placed FIRST in both panels (repeat --first for several, kept in the given order;
           the output folder is named after the first one) (e.g. the prostate WGBS models), given as a TSV with columns
           gene, cv_r2, probe (one row per selected CpG; "probe" = the 450K probe ID when the CpG is on the array, so it
           matches the TCGA probes, otherwise any other unique CpG ID such as chr_pos). Output goes to
           <WORK_DIR>/saturation_<file name>/. In the random-order panel it stays first and only the TCGA cohorts are permuted.
Writes <WORK_DIR>/saturation/saturation_plot.{pdf,png} and saturation_counts.tsv
"""
import sys, os, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta")
firsts = [a[len("--first="):] for a in sys.argv[2:] if a.startswith("--first=")]   # several allowed, kept in this order
first = firsts[0] if firsts else None
stable = next(((a.split("=", 1)[1] if "=" in a else "corr0.5") for a in sys.argv[2:] if a.startswith("--stable")), None)
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
CV_ALL = any(a == "--cv=all" for a in sys.argv[2:])   # --cv=all: no cv R2 filter, every model counts
TH = (-np.inf,) if CV_ALL else (0.1, 0.3); N_PERM = 500
thl = lambda th: "all models" if np.isinf(th) else f"cv R² ≥ {th}"
out = f"{work}/saturation" + ("_" + os.path.basename(first.rsplit("=", 1)[1].split(":")[0]).split(".")[0] if first else "") + next((a[len("--out_suffix="):] for a in sys.argv[2:] if a.startswith("--out_suffix=")), "") + ("_stable" if stable else "") + ("_cvall" if any(a == "--cv=all" for a in sys.argv[2:]) else ""); os.makedirs(out, exist_ok=True)

def stable_probes(c):
    """gene -> selected probes lying in a stable co-methylated region (as stability_report.py)."""
    S = f"{work}/{c}/stability_{SC}"
    cpg = pd.read_csv(f"{S}/cpg_freq.tsv", sep="\t", usecols=["TxID", "CpGID", "saved_selected"], dtype={"TxID": str})
    cpg = cpg[cpg.saved_selected == 1]
    rm = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 2, 3], names=["TxID", "CpGID", "region"], dtype={0: str})
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
        if not os.path.exists(sf): continue
        probes = [l.strip() for l in open(sf) if l.strip()]
        if st is not None:
            probes = [x for x in probes if x in st.get(g, ())]
            if not probes: continue
        rows.append((g, float(v["cv_r2"]), probes))
    return rows
data = {}
for c in cohorts:
    data[c] = load(c); print(f"{c}: {len(data[c]):,} gene models", flush=True)
nsamp = {}
for c in cohorts:
    f = f"{work}/{c}/need_gene_exp_{c}.csv"
    nsamp[c] = len(open(f).readline().strip().split(",")) - 1 if os.path.exists(f) else 0
n_tcga = len(cohorts)
FIRST = []
for fa in firsts:
    lab, f = fa.split("=", 1); n = 0
    if ":" in f: f, n = f.rsplit(":", 1)
    d = pd.read_csv(f, sep="\t")
    data[lab] = [(g, x.cv_r2.iloc[0], list(x.probe)) for g, x in d.groupby("gene")]; nsamp[lab] = int(n)
    FIRST.append(lab); print(f"{lab} (first): {len(data[lab]):,} gene models with >= 1 selected CpG", flush=True)
cohorts = FIRST + cohorts

# element ids: probes -> int; pairs -> gene_idx * N_PROBES + probe_idx. Each element gets a bitmask of the cohorts that
# contain it; the cumulative catalogue size for any cohort order is then computed from the (few) distinct masks.
P = {p: i for i, p in enumerate(sorted({x for c in cohorts for _, _, ps in data[c] for x in ps}))}
G = {g: i for i, g in enumerate(sorted({g for c in cohorts for g, _, _ in data[c]}))}
NP = len(P); ci = {c: k for k, c in enumerate(cohorts)}
def masks(kind, th):
    ids, bits = [], []
    for c in cohorts:
        e = [P[x] if kind == "probes" else G[g] * NP + P[x] for g, cv, ps in data[c] if cv >= th for x in ps]
        e = np.unique(np.array(e, dtype=np.int64)); ids.append(e); bits.append(np.full(len(e), 1 << ci[c], dtype=np.int64))
    ids, bits = np.concatenate(ids), np.concatenate(bits)
    if len(ids) == 0: return np.array([], dtype=np.int64), np.array([], dtype=np.int64)    # empty catalogue
    o = np.argsort(ids, kind="stable"); ids, bits = ids[o], bits[o]
    start = np.r_[0, np.flatnonzero(np.diff(ids)) + 1]
    m = np.bitwise_or.reduceat(bits, start)
    um, cnt = np.unique(m, return_counts=True)
    return um, cnt
def cumulative(um, cnt, order):
    rank = np.array([order.index(c) for c in cohorts]); first = np.full(len(um), len(cohorts))
    for k in range(len(cohorts)):
        first = np.where((um >> k) & 1, np.minimum(first, rank[k]), first)
    return np.array([cnt[first <= s].sum() for s in range(len(cohorts))]), np.array([cnt[first == s].sum() for s in range(len(cohorts))])

rng = np.random.default_rng(0); rows = []
U = {"probes": "CpGs", "pairs": "gene-CpG pairs"} if first else {"probes": "probes", "pairs": "pairs"}
fig, ax = plt.subplots(2, 2 * len(TH), figsize=(10 * len(TH), 9))
order = sorted(cohorts, key=lambda c: -nsamp[c])
if first: order = FIRST + [c for c in order if c not in FIRST]
rand_order = (lambda: FIRST + list(rng.permutation(cohorts[len(FIRST):]))) if first else (lambda: list(rng.permutation(cohorts)))
for j, (kind, th) in enumerate([(k, t) for k in ("probes", "pairs") for t in TH]):
    um, cnt = masks(kind, th)
    curves = np.array([cumulative(um, cnt, rand_order())[0] for _ in range(N_PERM)]); x = np.arange(1, len(cohorts) + 1)
    a = ax[0, j]
    a.fill_between(x, np.percentile(curves, 5, 0), np.percentile(curves, 95, 0), color="#1baf7a", alpha=.25, lw=0, label="5–95% of random orders" + (f" ({' + '.join(FIRST)} always first)" if first else ""))
    a.plot(x, curves.mean(0), "o-", color="#1f5f6b", ms=4, label="mean")
    a.set_xlabel("number of cohorts"); a.set_ylabel(f"distinct functional {U[kind]}")
    a.set_title(f"{U[kind]}, {thl(th)}: random cohort orders", loc="left", fontsize=10); a.legend(frameon=False, fontsize=8)
    cum, new = cumulative(um, cnt, order)
    own = {c: int(cnt[((um >> ci[c]) & 1).astype(bool)].sum()) for c in cohorts}
    a = ax[1, j]
    for i, c in enumerate(order):
        a.bar(i, cum[i] - new[i], color="#e1e0d9", width=.7); a.bar(i, new[i], bottom=cum[i] - new[i], color="#1baf7a", width=.7)
        rows.append(dict(catalogue=kind, cv_threshold=th, step=i + 1, cohort=c, n_samples=nsamp[c], cohort_size=own[c], new=int(new[i]),
                         cumulative=int(cum[i]), new_share_of_cohort=new[i] / max(1, own[c])))
    a.set_xticks(range(len(order)), [f"{c}\n({nsamp[c]})" for c in order], fontsize=7, rotation=90)
    a.set_ylabel(f"cumulative functional {U[kind]}"); a.set_title(f"{U[kind]}, {thl(th)}: by sample size (grey = already in, green = new)", loc="left", fontsize=10)
for a in ax.flat: a.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"Saturation of {'STABLE ' if stable else ''}functional CpGs across {(' + '.join(FIRST) + ' + ') if first else ''}{n_tcga} TCGA cohorts (METEOR, {SC} values; a CpG is functional if selected by a model{'' if CV_ALL else ' with cv R² ≥ threshold'})", fontsize=12)
plt.tight_layout(); plt.savefig(f"{out}/saturation_plot.pdf"); plt.savefig(f"{out}/saturation_plot.png", dpi=90)
R = pd.DataFrame(rows); R.to_csv(f"{out}/saturation_counts.tsv", sep="\t", index=False)
last = R[R.step == R.step.max()]
print(last[["catalogue", "cv_threshold", "cohort", "new", "new_share_of_cohort", "cumulative"]].to_string(index=False))
print(f"wrote {out}/saturation_plot.pdf, {out}/saturation_counts.tsv")
