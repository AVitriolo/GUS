"""Stable-only selected features for cluster-level stability selection (stability_selection_features.py +
stability_selection_summary.R): each model's saved selected features whose own selection frequency reaches the
gene's PFER cutoff (cpg_freq.tsv, stable = TRUE). Written in the sel_feats format (one file per model, lines
<TxID>_<feature>); models without any stable feature are not written.
Usage: python stable_sel_feats_features.py --stab_dir DIR --out_dir DIR --prefix <prefix>_<CV>
"""
import argparse
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--stab_dir", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--prefix", required=True)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

c = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", dtype={"CpGID": str})
s = c[(c["saved_selected"] == 1) & c["stable"]]
for tx, x in s.groupby("TxID"):
    with open(os.path.join(a.out_dir, f"{a.prefix}_{tx}"), "w") as fh:
        fh.write("".join(f"{tx}_{f}\n" for f in x["CpGID"]))
n_sel = int((c["saved_selected"] == 1).sum())
print(f"{s['TxID'].nunique():,} models with stable features; {len(s):,} of {n_sel:,} selected features stable -> {a.out_dir}")
