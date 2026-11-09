#!/usr/bin/env python
# =============================================================================
# test_cpg_localization.py
# Statistical test of the genomic localization of the selected CpGs
# (promoter / proximal / distal enhancer / exon / intron / intergenic, same
# exclusive categories and sources as annotate_selected_cpgs.py).
#
# A CpG-level Fisher test is invalid: neighbouring CpGs are correlated and are
# selected in clusters. Two tests that respect this:
#   1. per-model paired test: for each kept model (cv_r2 > mean), the fraction
#      of its selected CpGs in a category vs the fraction of its own candidate
#      CpGs (the +-1 Mb window); paired Wilcoxon signed-rank across models,
#      BH-adjusted over categories. The unit is the model, not the CpG.
#   2. inputs for a region-level permutation test (permtest_cpg_localization.R,
#      regioneR): one point per distinct selected region (selected CpGs of a
#      transcript chained <= gap bp, sel_feats_to_regions_w20.R), its middle
#      CpG; the universe is the candidate CpGs of the kept models.
#
# A model's candidates are the candidate CpGs between the first and last CpG of
# its window (identity map); the union of all windows shares one CpG filter.
# =============================================================================
import argparse, gzip, os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--features", required=True)      # TxID CpGID pos cv_r2 regionID coordID
ap.add_argument("--cvr2", required=True)          # TxID cv_r2 kept
ap.add_argument("--candidates", required=True)    # chr <TAB> pos (candidates of the kept transcripts)
ap.add_argument("--identity_dir", required=True)  # per-transcript candidate windows (*_<TxID>_identity_map.txt)
ap.add_argument("--identity_prefix", default="hg38_v38_pc_50000_25_10_")
ap.add_argument("--ccre", nargs="+", required=True)   # ENCODE cCRE beds, class in col 6 (PLS / pELS / dELS)
ap.add_argument("--gtf", required=True)
ap.add_argument("--out_dir", required=True)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
log = lambda *x: print(*x, flush=True)

# ---- annotation (as in annotate_selected_cpgs.py) ----
def merged(df):
    out = {}
    for ch, g in df.groupby("chr"):
        g = g.sort_values("start"); s = g.start.to_numpy(); e = g.end.to_numpy()
        new = np.r_[True, s[1:] > np.maximum.accumulate(e)[:-1]]
        out[ch] = (s[new], np.maximum.reduceat(e, np.flatnonzero(new)))
    return out
def hits(points, ivs):
    res = np.zeros(len(points), bool)
    for ch, g in points.groupby("chr"):
        if ch not in ivs: continue
        s, e = ivs[ch]; p = g.pos.to_numpy() - 1
        i = np.searchsorted(s, p, side="right") - 1
        res[g.index] = (i >= 0) & (p < e[np.clip(i, 0, None)])
    return res

cc = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chr", "start", "end", "cls"]) for f in a.ccre])
ann = {k: merged(cc[cc.cls == k]) for k in ["PLS", "pELS", "dELS"]}
gtf = pd.read_csv(a.gtf, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4], names=["chr", "feature", "start", "end"])
gtf["start"] -= 1
ann["exon"] = merged(gtf[gtf.feature == "exon"]); ann["gene"] = merged(gtf[gtf.feature == "gene"])
del cc, gtf
CATS = ["Promoter (PLS)", "Proximal enhancer (pELS)", "Distal enhancer (dELS)", "Exon", "Intron", "Intergenic"]
order = list(zip(["PLS", "pELS", "dELS", "exon", "intron", "intergenic"], CATS))

c = pd.read_csv(a.candidates, sep="\t", header=None, names=["chr", "pos"], dtype={"chr": str, "pos": np.int64})
c = c.drop_duplicates().sort_values(["chr", "pos"]).reset_index(drop=True)
h = {x: hits(c, ann[x]) for x in ["PLS", "pELS", "dELS", "exon", "gene"]}
h["intron"] = h["gene"] & ~h["exon"]; h["intergenic"] = ~h["gene"]
code = np.full(len(c), -1, np.int8)
for i, (x, _) in enumerate(order):
    code[(code == -1) & h[x]] = i
c["cat"] = code; del h
log(f"candidates: {len(c)} CpGs, unlabelled {int((code == -1).sum())}")

# ---- per-chromosome cumulative counts for window composition ----
chrom = {}
for ch, g in c.groupby("chr"):
    onehot = np.zeros((len(g) + 1, len(CATS)), np.int64)
    onehot[1:][np.arange(len(g)), g.cat.to_numpy()] = 1
    chrom[ch] = (g.pos.to_numpy(), np.cumsum(onehot, axis=0))
key = pd.Series(c.cat.to_numpy(), index=c.chr + "_" + c.pos.astype(str))

kept = pd.read_csv(a.cvr2, sep="\t").query("kept == True").TxID.tolist()
f = pd.read_csv(a.features, sep="\t", usecols=["TxID", "CpGID", "pos", "coordID"])
f = f[f.TxID.isin(set(kept))].copy()
f["cat"] = key.reindex(f.CpGID).to_numpy()
log(f"kept models: {len(kept)}, selected CpG rows: {len(f)}, not among candidates: {int(f.cat.isna().sum())}")

def first_last(path):
    with open(path, "rb") as fh:
        first = fh.readline().split(b"\t")[0].decode()
        fh.seek(-2, os.SEEK_END)
        while fh.read(1) != b"\n": fh.seek(-2, os.SEEK_CUR)
        last = fh.readline().split(b"\t")[0].decode()
    return first, last

sel_counts = f.dropna(subset=["cat"]).astype({"cat": int}).groupby(["TxID", "cat"]).size().unstack(fill_value=0)
sel_counts = sel_counts.reindex(columns=range(len(CATS)), fill_value=0)
rows = []
for tx in kept:
    lo, hi = first_last(os.path.join(a.identity_dir, f"{a.identity_prefix}{tx}_identity_map.txt"))
    ch = lo.rsplit("_", 1)[0]; p0 = int(lo.rsplit("_", 1)[1]); p1 = int(hi.rsplit("_", 1)[1])
    pos, cum = chrom[ch]
    i0 = np.searchsorted(pos, p0, side="left"); i1 = np.searchsorted(pos, p1, side="right")
    cand = cum[i1] - cum[i0]
    sel = sel_counts.loc[tx].to_numpy() if tx in sel_counts.index else np.zeros(len(CATS), int)
    rows.append([tx, *cand, *sel])
cols = [f"cand_{i}" for i in range(len(CATS))] + [f"sel_{i}" for i in range(len(CATS))]
m = pd.DataFrame(rows, columns=["TxID"] + cols)
m = m[(m.filter(like="sel_").sum(axis=1) > 0) & (m.filter(like="cand_").sum(axis=1) > 0)]
m.to_csv(os.path.join(a.out_dir, "per_model_category_counts.tsv.gz"), sep="\t", index=False, compression="gzip")

# ---- test 1: paired Wilcoxon across models ----
fs = m.filter(like="sel_").to_numpy() / m.filter(like="sel_").sum(axis=1).to_numpy()[:, None]
fc = m.filter(like="cand_").to_numpy() / m.filter(like="cand_").sum(axis=1).to_numpy()[:, None]
res = []
for i, lab in enumerate(CATS):
    d = fs[:, i] - fc[:, i]
    p = wilcoxon(fs[:, i], fc[:, i], zero_method="wilcox").pvalue
    pooled = fs[:, i].mean() / fc[:, i].mean()
    res.append([lab, len(d), 100 * fs[:, i].mean(), 100 * fc[:, i].mean(), pooled,
                100 * np.median(d), 100 * (d > 0).mean(), 100 * (d < 0).mean(), p])
r = pd.DataFrame(res, columns=["category", "n_models", "mean_pct_selected", "mean_pct_candidates",
                               "ratio_of_means", "median_diff_pct_points", "pct_models_higher",
                               "pct_models_lower", "p_wilcoxon"])
r["q_BH"] = multipletests(r.p_wilcoxon, method="fdr_bh")[1]
r.to_csv(os.path.join(a.out_dir, "per_model_wilcoxon.tsv"), sep="\t", index=False, float_format="%.4g")
log(r.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# ---- inputs for test 2: one point per distinct selected region ----
pts = (f.dropna(subset=["cat"]).sort_values(["coordID", "pos"])
        .groupby("coordID").apply(lambda g: g.iloc[(len(g) - 1) // 2][["CpGID", "pos", "cat"]]))
pts["chr"] = pts.CpGID.str.rsplit("_", n=1).str[0]
pts[["chr", "pos", "cat"]].astype({"cat": int}).to_csv(os.path.join(a.out_dir, "region_points.tsv.gz"), sep="\t", index=False, compression="gzip")
c[["chr", "pos", "cat"]].to_csv(os.path.join(a.out_dir, "candidate_categories.tsv.gz"), sep="\t", index=False, compression="gzip")
pd.Series(CATS).to_csv(os.path.join(a.out_dir, "categories.tsv"), sep="\t", index_label="cat", header=["category"])
log(f"region points: {len(pts)}; done")
