#!/usr/bin/env python
# =============================================================================
# stable_regions_annotation.py
# Genomic and repeat annotation of the CpGs in stable co-methylated regions
# (stability_selection_regions.R, corr regions), compared with the gene's other
# selected CpGs and with its candidate CpGs.
#   sets (per gene with >= 1 stable region):
#     stable     = CpGs selected by the model that lie in a stable region
#     unstable   = the model's other selected CpGs
#     candidates = candidate CpGs of the gene's 2Mb window (random n_cand per gene)
#   genomic : one category per CpG, first that applies (as annotate_selected_cpgs.py):
#             promoter (ENCODE PLS) > proximal enhancer (pELS) > distal enhancer
#             (dELS) > exon > intron > intergenic (GENCODE v38, all genes)
#   repeats : RepeatMasker class of the overlapping repeat with the lowest
#             divergence, "no repeat" otherwise
#     selected   = stable + unstable (all the model's selected CpGs)
# Test (as test_cpg_localization.py: the unit is the gene, not the CpG, since
# neighbouring CpGs are correlated and selected together): per gene, fraction of
# the set in the category; paired Wilcoxon signed-rank stable vs candidates,
# stable vs unstable and selected vs candidates across genes, BH over categories.
# A CpG is its C (0-based [pos - 1, pos)).
# =============================================================================
import argparse, glob, io, os, subprocess, tempfile
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--stab_dir", required=True)
ap.add_argument("--tag", default="corr0.5")
ap.add_argument("--cand_dir", required=True)
ap.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
ap.add_argument("--pls", required=True)
ap.add_argument("--els", required=True)
ap.add_argument("--gtf", required=True)
ap.add_argument("--rmsk", required=True)
ap.add_argument("--n_cand", type=int, default=1000, help="candidate CpGs sampled per gene")
ap.add_argument("--out_dir", required=True)
ap.add_argument("--reuse_table", action="store_true", help="reuse cpg_sets_annotation.tsv.gz, recompute statistics only")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
log = lambda *x: print(*x, flush=True)

CATS = ["Promoter (PLS)", "Proximal enhancer (pELS)", "Distal enhancer (dELS)", "Exon", "Intron", "Intergenic"]
table = os.path.join(a.out_dir, "cpg_sets_annotation.tsv.gz")
if a.reuse_table and os.path.exists(table):
    long = pd.read_csv(table, sep="\t")                     # annotation already done: statistics only
    tmp = None
else:
    # ---- per-gene sets ----
    reg = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t")
    stab = reg[reg["stable"]]
    genes = sorted(stab["TxID"].unique())
    cpg = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "saved_selected"])
    sel = cpg[(cpg["saved_selected"] == 1) & cpg["TxID"].isin(genes)][["TxID", "CpGID"]].copy()
    sel["pos"] = sel["CpGID"].str.split("_").str[1].astype(np.int64)
    sel = sel.merge(stab[["TxID", "start", "end"]], on="TxID", how="left")
    sel["in"] = (sel["pos"] >= sel["start"]) & (sel["pos"] <= sel["end"])
    sel = sel.groupby(["TxID", "CpGID"], as_index=False)["in"].any()
    sel["set"] = np.where(sel["in"], "stable", "unstable")
    def sample_cand(g):
        c = pd.read_csv(os.path.join(a.cand_dir, f"{a.prefix}_{g}.bed"), sep="\t", header=None,
                        names=["TxID", "CpGID"]).drop_duplicates()
        return c.sample(min(len(c), a.n_cand), random_state=1)
    cand = pd.concat([sample_cand(g) for g in genes])
    cand["set"] = "candidates"
    long = pd.concat([sel[["TxID", "CpGID", "set"]], cand], ignore_index=True)
    log(long.groupby("set").agg(pairs=("CpGID", "size"), distinct=("CpGID", "nunique")), f"\n{len(genes)} genes")

    # ---- annotate every distinct CpG once ----
    u = pd.DataFrame({"CpGID": long["CpGID"].unique()})
    u["chr"] = u["CpGID"].str.rsplit("_", n=1).str[0]
    u["pos"] = u["CpGID"].str.rsplit("_", n=1).str[1].astype(np.int64)


    def merged(df):
        out = {}
        for ch, g in df.groupby("chr"):
            g = g.sort_values("start"); s = g["start"].to_numpy(); e = g["end"].to_numpy()
            new = np.r_[True, s[1:] > np.maximum.accumulate(e)[:-1]]
            out[ch] = (s[new], np.maximum.reduceat(e, np.flatnonzero(new)))
        return out


    def hits(points, ivs):
        res = np.zeros(len(points), bool)
        for ch, g in points.groupby("chr"):
            if ch not in ivs:
                continue
            s, e = ivs[ch]; p = g["pos"].to_numpy() - 1
            i = np.searchsorted(s, p, side="right") - 1
            res[points.index.get_indexer(g.index)] = (i >= 0) & (p < e[np.clip(i, 0, None)])
        return res


    cc = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chr", "start", "end", "cls"])
                    for f in (a.pls, a.els)])
    gtf = pd.read_csv(a.gtf, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4], names=["chr", "feature", "start", "end"])
    gtf["start"] -= 1
    ann = {k: merged(cc[cc["cls"] == k]) for k in ["PLS", "pELS", "dELS"]}
    ann["exon"] = merged(gtf[gtf["feature"] == "exon"])
    ann["gene"] = merged(gtf[gtf["feature"] == "gene"])
    cat = np.full(len(u), "Intergenic", dtype=object)
    in_gene = hits(u, ann["gene"])
    cat[in_gene] = "Intron"
    for k, name in [("exon", "Exon"), ("dELS", "Distal enhancer (dELS)"), ("pELS", "Proximal enhancer (pELS)"),
                    ("PLS", "Promoter (PLS)")]:
        cat[hits(u, ann[k])] = name                                  # later = higher priority
    u["category"] = cat
    log("genomic categories done")

    tmp = tempfile.mkdtemp(dir=a.out_dir)
    rm = pd.read_csv(a.rmsk, sep="\t", header=None, usecols=[2, 5, 6, 7, 10, 11, 12],
                     names=["milliDiv", "chr", "start", "end", "repName", "repClass", "repFamily"])
    rm_bed, pt_bed = os.path.join(tmp, "rmsk.bed"), os.path.join(tmp, "points.bed")
    rm.sort_values(["chr", "start"])[["chr", "start", "end", "repClass", "milliDiv"]].to_csv(rm_bed, sep="\t", header=False, index=False)
    pd.DataFrame({"chr": u["chr"], "s": u["pos"] - 1, "e": u["pos"], "id": np.arange(len(u))}).sort_values(["chr", "s"]) \
        .to_csv(pt_bed, sep="\t", header=False, index=False)
    ov = subprocess.run(f"bedtools intersect -a {pt_bed} -b {rm_bed} -wa -wb -sorted", shell=True, check=True,
                        capture_output=True, text=True).stdout
    ov = pd.read_csv(io.StringIO(ov), sep="\t", header=None, usecols=[3, 7, 8], names=["id", "repClass", "milliDiv"])
    ov = ov.sort_values("milliDiv").drop_duplicates("id")
    u["repClass"] = "No repeat"
    u.loc[ov["id"].to_numpy(), "repClass"] = ov["repClass"].to_numpy()
    log("repeats done")

    long = long.merge(u[["CpGID", "category", "repClass"]], on="CpGID")
    long.to_csv(os.path.join(a.out_dir, "cpg_sets_annotation.tsv.gz"), sep="\t", index=False, compression="gzip")

# ---- pooled shares (distinct CpGs per set) and per-gene paired tests ----
# "selected" = all the model's selected CpGs (stable + unstable)
long = pd.concat([long[long["set"] != "selected"],
                  long[long["set"].isin(["stable", "unstable"])].assign(set="selected")], ignore_index=True)
SETS = ["candidates", "selected", "unstable", "stable"]
COMPARISONS = [("stable", "candidates"), ("stable", "unstable"), ("selected", "candidates")]
rows, tests = [], []
for col, cats in [("category", CATS), ("repClass", None)]:
    if cats is None:
        top = long[long["set"] == "candidates"][col].value_counts()
        cats = list(top.index[:8])
    for s in SETS:
        d = long[long["set"] == s].drop_duplicates("CpGID")
        for c in cats:
            rows.append({"annotation": col, "category": c, "set": s, "n": int((d[col] == c).sum()),
                         "pct": 100 * (d[col] == c).mean()})
    frac = long.groupby(["TxID", "set", col]).size().unstack(col, fill_value=0)
    frac = frac.div(frac.sum(axis=1), axis=0)
    for c in cats:
        f = frac[c].unstack("set") if c in frac else None
        if f is None:
            continue
        for s1, other in COMPARISONS:
            ff = f[[s1, other]].dropna()
            stat = wilcoxon(ff[s1], ff[other]) if (ff[s1] != ff[other]).any() else None
            tests.append({"annotation": col, "category": c, "comparison": f"{s1} vs {other}", "genes": len(ff),
                          "median_set": 100 * ff[s1].median(), "median_other": 100 * ff[other].median(),
                          "mean_diff_pct": 100 * (ff[s1] - ff[other]).mean(),
                          "p": stat.pvalue if stat else np.nan})
pooled = pd.DataFrame(rows)
tests = pd.DataFrame(tests)
tests["padj"] = np.nan
for key, g in tests.groupby(["annotation", "comparison"]):
    ok = g["p"].notna()
    tests.loc[g.index[ok], "padj"] = multipletests(g.loc[ok, "p"], method="fdr_bh")[1]
pooled.to_csv(os.path.join(a.out_dir, "annotation_pooled.tsv"), sep="\t", index=False, float_format="%.4g")
tests.to_csv(os.path.join(a.out_dir, "annotation_tests_per_gene.tsv"), sep="\t", index=False, float_format="%.4g")
if tmp:
    subprocess.run(["rm", "-rf", tmp])
log(pooled.pivot_table(index=["annotation", "category"], columns="set", values="pct").round(1)[SETS].to_string())
log(tests.round(4).to_string(index=False))
