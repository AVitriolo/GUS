"""CpG-level functional check of the MethylDriver CpG set with DAC in DU145 (GSE185651: EPIC GSE185594 + RNA-seq GSE185650).
EPIC: Ctl / DAC / DAC+triptolide (DACTPL), 5 d and 20 d, 2 reps (SWAN betas); RNA-seq: 20 d, 3 reps, authors' DESeq2.
Pairs: every candidate CpG of a modelled gene that is on EPIC (epic_candidate_index.tsv), train Pearson r with the gene from
test_replication_selected_cpgs.py --no_adjust; selected = sel_feats_..._stable_corr0.5.

Per gene, predicted response = sum over its selected EPIC CpGs of r_train * delta_beta (linearised: methylation loss at a
negatively associated CpG predicts a rise). Statistic: Spearman of that score with the gene's DESeq2 Wald stat across genes.
Nulls (n_draws): each selected CpG replaced by a non-selected candidate of the same gene with |r_train| within +-tol (nearest
if none) ("matched"), or by any non-selected candidate ("random"). Matched answers: do the selected CpGs predict better than
equally correlated CpGs of the same window? Secondary: % concordant sign(stat) == -sign(r) among demethylated pairs
(baseline beta >= 0.3, delta <= -0.1), selected vs matched draws.
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
p.add_argument("--pairs", required=True)
p.add_argument("--index", required=True)
p.add_argument("--betas", required=True)
p.add_argument("--rnaseq_dir", required=True)
p.add_argument("--rep_tx", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--n_draws", type=int, default=500)
p.add_argument("--tol", type=float, default=0.05)
a = p.parse_args()
rng = np.random.default_rng(1)

pr = pd.read_csv(a.pairs, sep="\t", usecols=["TxID", "CpGID", "set", "r_train"])
pr = pr.merge(pd.read_csv(a.index, sep="\t", usecols=["TxID", "CpGID", "Probe_ID"]), on=["TxID", "CpGID"])
b = pd.read_csv(a.betas, sep="\t", index_col=0)
mean = lambda c: b[[f"{c}_rep1", f"{c}_rep2"]].mean(1)
eb = pd.DataFrame({"beta_ctl": mean("Ctl_20d"), "db_DAC_20d": mean("DAC_20d") - mean("Ctl_20d"),
                   "db_DACTPL_20d": mean("DACTPL_20d") - mean("Ctl_20d"), "db_DAC_5d": mean("DAC_5d") - mean("Ctl_5d")})
pr = pr.merge(eb, left_on="Probe_ID", right_index=True)

rep = pd.read_csv(a.rep_tx, sep="\t")[["TxID", "gene"]]
rna = {}
for arm in ["DAC", "DACTPL"]:
    x = pd.read_csv(os.path.join(a.rnaseq_dir, f"rnaseq_DESeq2_{arm}_vs_Ctrl.tsv"), sep="\t")
    x.columns = ["refseq_id", "gene", "log2FC", "stat", "pvalue", "padj"]
    rna[arm] = rep.merge(x.dropna(subset=["stat"]).drop_duplicates("gene", keep=False), on="gene").set_index("TxID")
pr = pr[pr.TxID.isin(rna["DAC"].index)]
print(f"{pr.TxID.nunique():,} genes with RNA-seq; pairs: {pr.set.value_counts().to_dict()}")

# per gene: selected rows and, for each selected CpG, its |r|-matched control pool (indices into the gene's controls)
genes = []
for tx, g in pr.groupby("TxID"):
    s, c = g[g.set == "selected"], g[g.set == "random"]
    if len(s) == 0 or len(c) == 0:
        continue
    ca = c.r_train.abs().to_numpy()
    pools = []
    for rs in s.r_train.abs():
        w = np.flatnonzero(np.abs(ca - rs) <= a.tol)
        pools.append(w if len(w) else np.array([np.argmin(np.abs(ca - rs))]))
    genes.append((tx, s, c, pools))
print(f"{len(genes):,} genes with >= 1 selected and >= 1 control CpG on EPIC")


def score(r, db):
    return float(np.sum(r * db))


def run(db_col, arm, demeth=(0.3, -0.1)):
    stat = rna[arm].stat
    tx = [t for t, *_ in genes]
    y = stat.loc[tx].to_numpy()
    obs = np.array([score(s.r_train.to_numpy(), s[db_col].to_numpy()) for _, s, _, _ in genes])
    rho_obs = spearmanr(obs, y)[0]

    def conc(r, db, base, yy):
        m = (base >= demeth[0]) & (db <= demeth[1])
        return np.mean(np.sign(yy[m]) == -np.sign(r[m])) if m.any() else np.nan

    rep_y = np.concatenate([np.full(len(s), stat[t]) for t, s, _, _ in genes])
    cat = lambda f: np.concatenate([f(s) for _, s, _, _ in genes])
    conc_obs = conc(cat(lambda s: s.r_train.to_numpy()), cat(lambda s: s[db_col].to_numpy()),
                    cat(lambda s: s.beta_ctl.to_numpy()), rep_y)
    null = {"matched": [], "random": []}
    conc_null = []
    for _ in range(a.n_draws):
        for kind in null:
            sc, rr, dd, bb = [], [], [], []
            for _, s, c, pools in genes:
                if kind == "matched":
                    j = np.array([pl[rng.integers(len(pl))] for pl in pools])
                else:
                    j = rng.integers(len(c), size=len(s))
                cr, cd = c.r_train.to_numpy()[j], c[db_col].to_numpy()[j]
                sc.append(score(cr, cd))
                if kind == "matched":
                    rr.append(cr); dd.append(cd); bb.append(c.beta_ctl.to_numpy()[j])
            null[kind].append(spearmanr(sc, y)[0])
            if kind == "matched":
                conc_null.append(conc(np.concatenate(rr), np.concatenate(dd), np.concatenate(bb), rep_y))
    out = {"methylation": db_col, "rnaseq": f"{arm}_vs_Ctrl_20d", "n_genes": len(genes), "rho_selected": rho_obs,
           "p_spearman_selected": spearmanr(obs, y)[1]}
    for kind, v in null.items():
        v = np.array(v)
        out[f"rho_{kind}_mean"] = v.mean()
        out[f"rho_{kind}_95"] = np.quantile(v, 0.95)
        out[f"p_vs_{kind}"] = (1 + np.sum(v >= rho_obs)) / (1 + len(v))
    cn = np.array(conc_null)
    out.update({"pct_concordant_selected": 100 * conc_obs, "pct_concordant_matched_mean": 100 * np.nanmean(cn),
                "p_concordance_vs_matched": (1 + np.sum(cn >= conc_obs)) / (1 + len(cn))})
    return out, null, obs, y


res, keep = [], {}
for db_col, arm in [("db_DAC_20d", "DAC"), ("db_DACTPL_20d", "DACTPL"), ("db_DAC_5d", "DAC")]:
    o, null, obs, y = run(db_col, arm)
    res.append(o)
    keep[(db_col, arm)] = (o, null, obs, y)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items()}, flush=True)
res = pd.DataFrame(res)
res.to_csv(os.path.join(a.out_dir, "cpg_level_dac_du145.tsv"), sep="\t", index=False)

o, null, obs, y = keep[("db_DAC_20d", "DAC")]
pd.DataFrame({"TxID": [t for t, *_ in genes], "n_selected_epic": [len(s) for _, s, _, _ in genes],
              "pred_score": obs, "deseq2_stat": y}).to_csv(os.path.join(a.out_dir, "gene_scores_DAC_20d.tsv"), sep="\t", index=False)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
for kind, colr in [("random", "#bababa"), ("matched", "#67a9cf")]:
    ax[0].hist(null[kind], bins=30, color=colr, alpha=0.9, label=f"{kind} non-selected CpGs")
ax[0].axvline(o["rho_selected"], color="#b2182b", lw=2, label="selected (MethylDriver set)")
ax[0].set_xlabel("Spearman(predicted score, DESeq2 stat) across genes")
ax[0].set_ylabel(f"draws (n={a.n_draws})")
ax[0].legend(fontsize=7, frameon=False)
ax[0].set_title("DU145, DAC 20 d: EPIC delta beta vs RNA-seq", fontsize=9)
q = pd.qcut(pd.Series(obs).rank(method="first"), 10, labels=False)
m = pd.Series(y).groupby(q).median()
ax[1].plot(m.index + 1, m.values, "o-", color="#b2182b")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set_xlabel("decile of predicted score  sum(r_train x delta beta)")
ax[1].set_ylabel("median DESeq2 Wald stat (DAC vs Ctrl)")
ax[1].set_title("genes ranked by predicted response", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(a.out_dir, "cpg_level_dac_du145.png"), dpi=200)
fig.savefig(os.path.join(a.out_dir, "cpg_level_dac_du145.pdf"))
