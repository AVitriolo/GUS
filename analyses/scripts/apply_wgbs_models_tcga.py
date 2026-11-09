#!/usr/bin/env python
# Scores the EXACT WGBS-trained pipelines (scaler -> SelectKBest -> SVR) from
# run_tcga_transfer_hm450_2mb_kcap_vs_uncapped.sh on TCGA PRAD (HM450 + RNA-seq),
# uncapped vs k-capped, to see which one generalizes out of cohort. Nothing is refit.
#
# Feature mapping: WGBS features are "<TxID>_chrN_<pos>" (hg38); pos == HM450 manifest
# CpG_beg + 1, so each maps to one probe ID. PRAD_data_meth.csv is ALREADY M-values
# (read.data.preprocessing_function_infinite.r: log2(b/(1-b)), +-10 at b = 1/0, probes
# >20% NA dropped, the rest imputed), the same scale as the WGBS M assay - used as is.
# Training candidates were restricted to probes present there, so a missing feature is
# an error; any residual NA is filled with the TCGA per-probe median (label-free).
#
# Target: TxID-level kallisto (WGBS cohort) vs gene-level TCGA expression, on different
# scales - so rank/linear agreement (Spearman, Pearson) is the primary metric; raw R2
# (uncalibrated, scale shift included) is reported too.
import argparse
import os
import subprocess
from io import StringIO

import joblib as jbl
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import r2_score

p = argparse.ArgumentParser()
p.add_argument("--txid_list", required=True, help="TSV with TxID and gene columns")
p.add_argument("--results_dir", required=True, help="dir with {uncapped,kcapped}/{models,hyper,performance}")
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
p.add_argument("--cv", default=3, type=int)
p.add_argument("--manifest", default="resources/illumina_probes/HM450.hg38.manifest.tsv.gz")
p.add_argument("--tcga_meth", default="/home/ieo7293/TCGA_PRAD/PRAD_data_meth.csv")
p.add_argument("--tcga_expr", default="/home/ieo7293/TCGA_PRAD/need_gene_exp_PRAD.csv")
p.add_argument("--out_per_model", required=True)
p.add_argument("--out_summary", required=True)
args = p.parse_args()

CONDS = ["uncapped", "kcapped"]

txids = pd.read_csv(args.txid_list, sep="\t")

man = pd.read_csv(args.manifest, sep="\t", usecols=["CpG_chrm", "CpG_beg", "Probe_ID"]).dropna()
man["key"] = man["CpG_chrm"] + "_" + (man["CpG_beg"].astype(int) + 1).astype(str)
key2probe = dict(zip(man["key"], man["Probe_ID"]))


def strip_tx(feat, tx):
    return feat[len(tx) + 1:] if feat.startswith(tx + "_") else feat


# ---- load models + feature lists ----
models, feat_probe = {}, {}
for tx in txids.TxID:
    for c in CONDS:
        pkl = [f for f in os.listdir(f"{args.results_dir}/{c}/models") if f.startswith(f"{tx}_") and f.endswith(".pkl")]
        assert len(pkl) == 1, f"{tx} {c}: expected 1 pkl, found {pkl}"
        models[(tx, c)] = jbl.load(f"{args.results_dir}/{c}/models/{pkl[0]}")
    feats = list(models[(tx, CONDS[0])].feature_names_in_)
    assert feats == list(models[(tx, CONDS[1])].feature_names_in_), f"{tx}: feature order differs between conditions"
    feat_probe[tx] = {f: key2probe.get(strip_tx(f, tx)) for f in feats}

# ---- one pass over the 3.7 GB TCGA methylation CSV for all needed probes ----
probes = sorted({pr for d in feat_probe.values() for pr in d.values() if pr})
header = pd.read_csv(args.tcga_meth, nrows=0).columns.tolist()
pat = "\n".join(f'"{pr}",' for pr in probes)
hits = subprocess.run(["grep", "-F", "-f", "/dev/stdin", args.tcga_meth], input=pat,
                      capture_output=True, text=True, check=False).stdout
meth = pd.read_csv(StringIO(hits), header=None, names=header, index_col=0).astype("float64")
meth = meth.loc[meth.index.intersection(probes)]
meth_M = meth.apply(lambda r: r.fillna(r.median()), axis=1)
print(f"TCGA: {len(meth_M)}/{len(probes)} needed probes found")

expr = pd.read_csv(args.tcga_expr, index_col=0)
samples = sorted(set(meth_M.columns) & set(expr.columns))
tumour = [s for s in samples if s[13:15] == "01"]
print(f"TCGA samples with meth+expr: {len(samples)} ({len(tumour)} primary tumour)")


def read_kv(path):
    return pd.read_csv(path, sep="\t", header=None, index_col=0)[1] if os.path.exists(path) else pd.Series(dtype=object)


def metrics(y, yhat, tag):
    if np.std(yhat) == 0:
        return {f"{tag}_spearman": np.nan, f"{tag}_pearson": np.nan, f"{tag}_r2_raw": r2_score(y, yhat)}
    return {f"{tag}_spearman": stats.spearmanr(y, yhat)[0],
            f"{tag}_pearson": stats.pearsonr(y, yhat)[0],
            f"{tag}_r2_raw": r2_score(y, yhat)}


rows = []
for tx, gene in zip(txids.TxID, txids.gene):
    if gene not in expr.index:
        print(f"[{tx}] {gene} not in TCGA expression - skipped")
        continue
    fp = feat_probe[tx]
    missing = [f for f, pr in fp.items() if pr is None or pr not in meth_M.index]
    assert not missing, f"[{tx}] {len(missing)} features not measured in TCGA, e.g. {missing[:3]}"
    X = pd.DataFrame({f: meth_M.loc[pr, samples].values for f, pr in fp.items()}, index=samples)
    y = expr.loc[gene, samples].astype(float)

    for c in CONDS:
        pipe = models[(tx, c)]
        sel = set(np.array(pipe.feature_names_in_)[pipe.named_steps["feature_selection"].get_support()])
        yhat = pipe.predict(X)
        perf = read_kv(f"{args.results_dir}/{c}/performance/{args.prefix}_{args.cv}_{tx}")
        hyp = read_kv(f"{args.results_dir}/{c}/hyper/{args.prefix}_{args.cv}_{tx}")
        row = {"TxID": tx, "gene": gene, "condition": c, "n_feat": len(fp),
               "k": int(float(hyp.get("feature_selection__k", np.nan))),
               "kernel": hyp.get("model__kernel"),
               "n_selected": len(sel),
               "wgbs_cv_r2": float(perf.get("cv_r2", np.nan)),
               "wgbs_test_r2": float(perf.get("test_r2", np.nan))}
        row.update(metrics(y.values, yhat, "tcga_all"))
        ti = [samples.index(s) for s in tumour]
        row.update(metrics(y.values[ti], yhat[ti], "tcga_tumour"))
        rows.append(row)

res = pd.DataFrame(rows)
res.to_csv(args.out_per_model, sep="\t", index=False)

# ---- paired uncapped vs kcapped summary ----
summ = []
wide_cols = ["k", "wgbs_cv_r2", "wgbs_test_r2", "tcga_all_spearman", "tcga_all_pearson",
             "tcga_all_r2_raw", "tcga_tumour_spearman", "tcga_tumour_pearson"]
w = res.pivot(index="TxID", columns="condition", values=wide_cols)
for m in wide_cols:
    a, b = w[(m, "uncapped")].astype(float), w[(m, "kcapped")].astype(float)
    ok = a.notna() & b.notna()
    d = (b - a)[ok]
    pval = stats.wilcoxon(b[ok], a[ok]).pvalue if (d != 0).sum() > 0 else np.nan
    summ.append({"metric": m, "n": int(ok.sum()),
                 "median_uncapped": a[ok].median(), "median_kcapped": b[ok].median(),
                 "median_diff_kcapped_minus_uncapped": d.median(),
                 "n_kcapped_better": int((d > 0).sum()), "n_uncapped_better": int((d < 0).sum()),
                 "wilcoxon_p": pval})
summ = pd.DataFrame(summ)
summ.to_csv(args.out_summary, sep="\t", index=False)
print(summ.to_string(index=False))
