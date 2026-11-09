"""Export the stable functional-CpG catalogue of the TCGA cohorts as a table (one row per gene model x probe x cohort).
Same selection as saturation_plot.py --stable --cv=all: a probe counts for a gene model when the model selected it and it
lies in a stable co-methylated region (stability selection, region_freq_<TAG>.tsv). No cv R2 filter; cv_r2 is kept.
Usage: python scripts/export_stable_catalogue.py <WORK_DIR> [--stable=corr0.5] [--out=FILE.tsv.gz] [COHORT ...]
       (default: every cohort with results; METH_SCALE as in saturation_plot.py, default beta)
Output columns: cohort, gene, cv_r2, probe
"""
import sys, os, glob
import pandas as pd

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta")
stable = next((a.split("=", 1)[1] for a in sys.argv[2:] if a.startswith("--stable=")), "corr0.5")
out = next((a.split("=", 1)[1] for a in sys.argv[2:] if a.startswith("--out=")), f"{work}/stable_catalogue_{stable}_cvall.tsv.gz")
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(
    os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))

def stable_probes(c):
    """gene -> selected probes lying in a stable co-methylated region (as saturation_plot.py / stability_report.py)."""
    S = f"{work}/{c}/stability_{SC}"
    cpg = pd.read_csv(f"{S}/cpg_freq.tsv", sep="\t", usecols=["TxID", "CpGID", "saved_selected"], dtype={"TxID": str})
    cpg = cpg[cpg.saved_selected == 1]
    rm = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 2, 3], names=["TxID", "CpGID", "region"], dtype={0: str})
                    for f in glob.glob(f"{S}/region_map_{stable}_shard*.tsv")]).drop_duplicates(["TxID", "CpGID"])
    rf = pd.read_csv(f"{S}/region_freq_{stable}.tsv", sep="\t", usecols=["TxID", "region", "stable"], dtype={"TxID": str})
    m = cpg.merge(rm, on=["TxID", "CpGID"]).merge(rf[rf.stable], on=["TxID", "region"])
    return m.groupby("TxID").CpGID.apply(set).to_dict()

rows = []
for c in cohorts:
    st = stable_probes(c); n = 0
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        g = os.path.basename(f)[:-8]
        v = dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())
        sf = f"{work}/{c}/results_{SC}/sel_feats/{g}_SVR.tsv"
        if not os.path.exists(sf): continue
        probes = [x for x in (l.strip() for l in open(sf)) if x and x in st.get(g, ())]
        rows += [(c, g, float(v["cv_r2"]), p) for p in probes]; n += bool(probes)
    print(f"{c}: {n:,} gene models with >= 1 stable probe", flush=True)

d = pd.DataFrame(rows, columns=["cohort", "gene", "cv_r2", "probe"])
d.to_csv(out, sep="\t", index=False)
print(f"{len(d):,} rows, {d.probe.nunique():,} distinct probes, {d.cohort.nunique()} cohorts -> {out}")
