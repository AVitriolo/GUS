"""CpG-specific check of the MethylDriver CpG set with DAC WGBS + RNA-seq (GSE331134, part of GSE331137): MDS-L (DMSO,
0.5/1/2 uM), MOLM13 (DMSO, 1 uM), MOLM14 (DMSO, 1/2.5/3.5/5 uM); one sample per condition, per-CpG ratios without coverage.

DAC demethylates nearly every CpG, so "the CpG lost methylation" cannot single out the selected CpGs. What varies is HOW MUCH
each CpG lost. Per gene g and treated condition c:
    loss = (beta_DMSO - beta_DAC) / beta_DMSO, over CpGs with beta_DMSO >= --min_beta in that line
    L_sel  = mean loss of g's selected CpGs;  L_bg = mean loss of g's other candidate CpGs (its 2Mb window)
    excess = L_sel - L_bg   (> 0: the selected CpGs lost more than their window; < 0: they resisted)
    y      = s_g * log2FC(CPM) DAC vs DMSO, s_g = -sign(mean train r of the selected CpGs) (y > 0 = predicted direction)
Statistic: Spearman(excess, y) over gene-condition points (conditions pooled, each centred). If the selected CpGs matter
beyond their window, it is > 0. Null: the same with, per gene, an equal-size random set of non-selected candidates in place
of the selected CpGs (n_draws). Also y by excess tertile (selected demethylated more / as much / less than the window).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--wgbs_dir", required=True)
p.add_argument("--counts", required=True)
p.add_argument("--cand", required=True, help="cpg_freq.tsv (TxID, CpGID): candidate CpGs per gene")
p.add_argument("--sel_r", required=True, help="selected pairs with r_train (TxID, CpGID, set, r_train)")
p.add_argument("--sel_pairs", required=True)
p.add_argument("--rep_tx", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--min_beta", type=float, default=0.3)
p.add_argument("--min_cpm", type=float, default=1.0)
p.add_argument("--n_draws", type=int, default=200)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(1)

# conditions: WGBS file stem -> (line, dose label) and the matching RNA-seq column
COND = [("MDSL_DMSO_5d", "MDS-L", "0", "MDS-L_0uM"), ("MDSL_DAC0.5uM_5d", "MDS-L", "0.5", "MDS-L_0.5uM"),
        ("MDSL_DAC1.0uM_5d", "MDS-L", "1", "MDS-L_1uM"), ("MDSL_DAC2.0uM_5d", "MDS-L", "2", "MDS-L_2uM"),
        ("MOLM13_DMSO_5d", "MOLM13", "0", "MOLM13_0uM"), ("MOLM13_DAC1.0uM_5d", "MOLM13", "1", "MOLM13_1uM"),
        ("MOLM14_DMSO_5d", "MOLM14", "0", "MOLM14_0uM"), ("MOLM14_DAC1.0uM_5d", "MOLM14", "1", "MOLM14_1uM"),
        ("MOLM14_DAC2.5uM_5d", "MOLM14", "2.5", "MOLM14_2.5uM"), ("MOLM14_DAC3.5uM_5d", "MOLM14", "3.5", "MOLM14_3.5uM"),
        ("MOLM14_DAC5.0uM_5d", "MOLM14", "5", "MOLM14_5uM_5d")]

# ---- candidate pairs (gene x CpG), selected flag, CpG integer index ----
cand = pd.read_csv(a.cand, sep="\t", usecols=["TxID", "CpGID"])
sel = pd.read_csv(a.sel_pairs, sep="\t")[["TxID", "CpGID"]].assign(sel=True)
cand = cand.merge(sel, on=["TxID", "CpGID"], how="left").fillna({"sel": False})
cpg = pd.Index(cand.CpGID.unique())
cand["ci"] = cpg.get_indexer(cand.CpGID)
cand["gi"], genes_tx = pd.factorize(cand.TxID)
print(f"{len(cand):,} candidate pairs, {len(cpg):,} CpGs, {len(genes_tx):,} genes, {int(cand.sel.sum()):,} selected", flush=True)

# ---- WGBS: their chr:pos is our chr_pos + 1 ----
B = np.full((len(COND), len(cpg)), np.nan, dtype=np.float32)
for k, (stem, *_ ) in enumerate(COND):
    f = [x for x in os.listdir(a.wgbs_dir) if x.endswith(f"{stem}.CpG_meth_ratio.txt.gz")][0]
    x = pd.read_csv(os.path.join(a.wgbs_dir, f), sep="\t", header=None, names=["k", "v"])
    cp = x.k.str.split(":", n=1, expand=True)
    idx = cpg.get_indexer(cp[0] + "_" + (cp[1].astype(int) - 1).astype(str))
    m = idx >= 0
    B[k, idx[m]] = x.v.to_numpy()[m]
    print(f"{stem}: {m.sum():,} candidate CpGs covered", flush=True)

# ---- RNA: exon rows summed per gene name, CPM ----
cn = pd.read_csv(a.counts, sep="\t")
cn = cn.groupby("gene_name")[[c[3] for c in COND]].sum()
cpm = cn / cn.sum() * 1e6
rep = pd.read_csv(a.rep_tx, sep="\t")[["TxID", "gene"]].set_index("TxID").gene

# ---- per-gene direction from the selected CpGs ----
r = pd.read_csv(a.sel_r, sep="\t", usecols=["TxID", "CpGID", "set", "r_train"])
r = r[r.set == "selected"].groupby("TxID").r_train.agg(mean_r="mean", frac_neg=lambda v: (v < 0).mean())
s_g = -np.sign(r.mean_r.reindex(genes_tx)).to_numpy()
clear = ((r.frac_neg >= 0.75) | (r.frac_neg <= 0.25)).reindex(genes_tx).fillna(False).to_numpy()

gi, ci, is_sel = cand.gi.to_numpy(), cand.ci.to_numpy(), cand.sel.to_numpy()
G = len(genes_tx)
nsel_g = np.bincount(gi[is_sel], minlength=G)


def gene_mean(loss_pair, mask):
    ok = mask & np.isfinite(loss_pair)
    n = np.bincount(gi[ok], minlength=G)
    s = np.bincount(gi[ok], weights=loss_pair[ok], minlength=G)
    return np.where(n >= 3, s / np.maximum(n, 1), np.nan)


rows, pts = [], []
lines = {}
for k, (stem, line, dose, rcol) in enumerate(COND):
    lines.setdefault(line, []).append(k)
treated = [(ks[0], k, line, COND[k][2], COND[k][3]) for line, ks in lines.items() for k in ks[1:]]

loss_by_cond = {}
for k0, k, line, dose, rcol in treated:
    b0, b1 = B[k0][ci], B[k][ci]
    loss = np.where(b0 >= a.min_beta, (b0 - b1) / b0, np.nan)
    loss_by_cond[k] = loss
    c0 = cpm[COND[k0][3]]
    expr = (c0 >= a.min_cpm) | (cpm[rcol] >= a.min_cpm)
    lfc = np.log2(cpm[rcol] + 1) - np.log2(c0 + 1)
    gname = rep.reindex(genes_tx).to_numpy()
    y = pd.Series(lfc[expr]).reindex(gname).to_numpy() * s_g
    pts.append(pd.DataFrame({"TxID": genes_tx, "line": line, "dose": dose, "k": k, "y": y,
                             "L_sel": gene_mean(loss, is_sel), "L_bg": gene_mean(loss, ~is_sel), "clear": clear}))
pts = pd.concat(pts, ignore_index=True)
pts["excess"] = pts.L_sel - pts.L_bg
ok = pts.dropna(subset=["y", "excess"])
print(f"{len(ok):,} gene-condition points ({ok.TxID.nunique():,} genes)", flush=True)


def stat(d, col="excess"):
    d = d.dropna(subset=["y", col])
    cen = lambda c: d[c] - d.groupby("k")[c].transform("mean")
    return spearmanr(cen(col), cen("y"))[0]


obs = {"all genes": stat(ok), "clear-direction genes": stat(ok[ok.clear])}
print("observed:", {k: round(v, 4) for k, v in obs.items()}, flush=True)

# ---- null: random non-selected sets of the gene's selected size ----
ns_idx = np.flatnonzero(~is_sel)
null = {k: [] for k in obs}
for d in range(a.n_draws):
    key = rng.random(len(ns_idx))
    order = np.lexsort((key, gi[ns_idx]))
    srt = ns_idx[order]
    gs = gi[srt]
    start = np.searchsorted(gs, np.arange(G))
    rank = np.arange(len(srt)) - start[gs]
    pick = np.zeros(len(gi), bool)
    pick[srt[rank < nsel_g[gs]]] = True
    pz = []
    for k0, k, line, dose, rcol in treated:
        loss = loss_by_cond[k]
        pz.append(gene_mean(loss, pick) - gene_mean(loss, ~is_sel & ~pick))
    q = pts.copy()
    q["excess_null"] = np.concatenate(pz)
    q = q.dropna(subset=["y", "excess_null"])
    null["all genes"].append(stat(q, "excess_null"))
    null["clear-direction genes"].append(stat(q[q.clear], "excess_null"))
    if (d + 1) % 25 == 0:
        print(f"draw {d + 1}/{a.n_draws}", flush=True)

res = []
for kname, v in obs.items():
    nv = np.array(null[kname])
    res.append({"genes": kname, "rho_selected": v, "rho_null_mean": nv.mean(), "rho_null_sd": nv.std(),
                "rho_null_95": np.quantile(nv, 0.95), "p_empirical": (1 + np.sum(nv >= v)) / (1 + len(nv))})
res = pd.DataFrame(res)

# per condition and tertile summaries (clear-direction genes)
cc = ok[ok.clear].copy()
cc["tertile"] = cc.groupby("k").excess.transform(lambda e: pd.qcut(e, 3, labels=["resisted vs window", "as window", "lost more than window"]))
ter = cc.groupby("tertile", observed=True).agg(points=("y", "size"), mean_y=("y", "mean"), median_y=("y", "median"),
                                              pct_predicted_dir_gt0_5=("y", lambda v: 100 * np.mean(v > 0.5)),
                                              pct_opposite_gt0_5=("y", lambda v: 100 * np.mean(v < -0.5)))
percond = ok.groupby(["line", "dose"]).apply(lambda d: pd.Series({"points": len(d), "median_L_sel": d.L_sel.median(),
          "median_L_bg": d.L_bg.median(), "rho": spearmanr(d.excess, d.y)[0]}))
pts.to_csv(os.path.join(a.out_dir, "gene_condition_points.tsv"), sep="\t", index=False)
res.to_csv(os.path.join(a.out_dir, "cpg_specific_test.tsv"), sep="\t", index=False)
ter.to_csv(os.path.join(a.out_dir, "tertiles_clear_genes.tsv"), sep="\t")
percond.to_csv(os.path.join(a.out_dir, "per_condition.tsv"), sep="\t")
with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
    print(res, ter, percond, sep="\n\n")

fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ax[0].hist(null["clear-direction genes"], bins=30, color="#bababa", label="random non-selected sets")
ax[0].axvline(obs["clear-direction genes"], color="#b2182b", lw=2, label="selected (MethylDriver set)")
ax[0].set_xlabel("Spearman(excess demethylation, response in predicted direction)")
ax[0].set_ylabel(f"draws (n={a.n_draws})")
ax[0].legend(fontsize=7, frameon=False)
ax[0].set_title("AML lines, DAC: selected vs window", fontsize=9)
ax[1].bar(range(len(ter)), ter.mean_y, color=["#67a9cf", "#bababa", "#b2182b"])
ax[1].set_xticks(range(len(ter)), [t.replace(" vs window", "\nvs window").replace(" than window", "\nthan window") for t in ter.index], fontsize=8)
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set_ylabel("mean log2FC in predicted direction")
ax[1].set_title("selected CpGs' demethylation relative to window", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(a.out_dir, "cpg_specific_test.png"), dpi=200)
fig.savefig(os.path.join(a.out_dir, "cpg_specific_test.pdf"))
