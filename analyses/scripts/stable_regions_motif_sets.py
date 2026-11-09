#!/usr/bin/env python
# =============================================================================
# stable_regions_motif_sets.py
# Foreground / background BEDs for known-motif enrichment of the stable
# co-methylated regions (stability_selection_regions.R, corr regions).
#   foreground : stable regions holding >= 1 of the model's selected CpGs,
#                all genes pooled, overlapping regions merged; first - last CpG,
#                or (--fixed_len) a fixed window around the region centre;
#                padded around the centre to >= min_len bp
#   background : candidate co-methylated regions of the same gene windows
#                (bg_regions_corr<thr>_shard*.tsv, sampled regardless of
#                selection), not overlapping the foreground, padded the same way,
#                matched to the foreground on length, GC and CpG density
#                (strata = bins of the three), ratio x foreground per stratum
# Output BEDs use Ensembl chromosome names (GRCh38.p14 FASTA), as
# gimme_scan.sh; background_match.tsv summarises the matching.
# =============================================================================
import argparse, glob, os
import numpy as np
import pandas as pd
from pyfaidx import Fasta

ap = argparse.ArgumentParser()
ap.add_argument("--stab_dir", required=True)
ap.add_argument("--tag", default="corr0.5")
ap.add_argument("--genome_fa", required=True)
ap.add_argument("--min_len", type=int, default=100)
ap.add_argument("--fixed_len", type=int, default=0,
                help="if > 0, every region becomes a window of this size around its centre (equal lengths: "
                     "best-match motif scores grow with sequence length)")
ap.add_argument("--ratio", type=int, default=2)
ap.add_argument("--out_dir", required=True)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(1)
fa = Fasta(a.genome_fa, as_raw=True, sequence_always_upper=True)


def merge(df):
    """Merge overlapping / touching intervals (BED, 0-based half-open) per chromosome."""
    out = []
    for ch, g in df.sort_values(["chr", "start"]).groupby("chr", sort=False):
        s, e = g["start"].to_numpy(), g["end"].to_numpy()
        e_run = np.maximum.accumulate(e)
        new = np.r_[True, s[1:] > e_run[:-1]]
        out.append(pd.DataFrame({"chr": ch, "start": s[new], "end": np.maximum.reduceat(e, np.flatnonzero(new))}))
    return pd.concat(out, ignore_index=True)


def pad(df):
    if a.fixed_len:
        mid = (df["start"] + df["end"]) // 2
        df["start"] = (mid - a.fixed_len // 2).clip(lower=0)
        df["end"] = df["start"] + a.fixed_len
        return df
    ln = df["end"] - df["start"]
    short = ln < a.min_len
    mid = (df["start"] + df["end"]) // 2
    df.loc[short, "start"] = (mid - a.min_len // 2)[short]
    df.loc[short, "end"] = (df["start"] + a.min_len)[short]
    df["start"] = df["start"].clip(lower=0)
    return df


def features(df):
    gc, cpg = [], []
    for ch, s, e in zip(df["chr"], df["start"], df["end"]):
        seq = fa[ch.removeprefix("chr")][s:e]
        n = max(len(seq) - seq.count("N"), 1)
        gc.append((seq.count("G") + seq.count("C")) / n)
        cpg.append(seq.count("CG") / n * 100)
    df["len"] = df["end"] - df["start"]
    df["gc"], df["cpg_per100"] = gc, cpg
    return df


# ---- foreground: stable regions with >= 1 selected CpG (CpG positions 1-based -> BED) ----
r = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t")
fg = r[r["stable"] & (r["n_saved"] > 0)]
fg = merge(pd.DataFrame({"chr": fg["chr"], "start": fg["start"] - 1, "end": fg["end"]}))
fg = features(pad(fg))
print(f"foreground: {len(fg)} merged stable regions")

# ---- background pool: candidate regions of the same windows, not overlapping the foreground ----
bg = pd.concat([pd.read_csv(f, sep="\t", header=None, names=["TxID", "chr", "start", "end", "n_cpg"])
                for f in glob.glob(os.path.join(a.stab_dir, f"bg_regions_{a.tag}_shard*.tsv"))])
bg = pad(pd.DataFrame({"chr": bg["chr"], "start": bg["start"] - 1, "end": bg["end"]}).drop_duplicates()
         .reset_index(drop=True))
fgm = merge(fg[["chr", "start", "end"]].copy())
keep = np.ones(len(bg), bool)
for ch, idx in bg.groupby("chr").groups.items():
    f = fgm[fgm["chr"] == ch]
    if f.empty:
        continue
    s, e = f["start"].to_numpy(), f["end"].to_numpy()
    b = bg.loc[idx]
    i = np.searchsorted(s, b["end"].to_numpy(), side="left") - 1          # last fg starting before bg end
    keep[idx] =~((i >= 0) & (e[np.clip(i, 0, None)] > b["start"].to_numpy()))
bg = features(bg[keep].reset_index(drop=True))
print(f"background pool: {len(bg)} candidate regions not overlapping the foreground")

# ---- match on length, GC, CpG density ----
len_edges = np.unique(np.quantile(fg["len"], np.linspace(0, 1, 11)))
gc_edges = np.arange(0, 1.0001, 0.05)
cpg_edges = np.unique(np.r_[0, np.quantile(fg["cpg_per100"], np.linspace(0.1, 1, 10))])
for d in (fg, bg):
    d["stratum"] = (pd.Series(np.digitize(d["len"], len_edges[1:-1]), index=d.index).astype(str) + "_" +
                    pd.Series(np.digitize(d["gc"], gc_edges[1:-1]), index=d.index).astype(str) + "_" +
                    pd.Series(np.digitize(d["cpg_per100"], cpg_edges[1:-1]), index=d.index).astype(str))
pick, short = [], 0
pool = bg.groupby("stratum").indices
for st, n in fg["stratum"].value_counts().items():
    cand = pool.get(st, np.array([], int))
    want = a.ratio * n
    short += max(0, want - len(cand))
    pick.append(rng.choice(cand, size=min(want, len(cand)), replace=False))
bgm = bg.loc[np.concatenate(pick)]
print(f"matched background: {len(bgm)} regions (short by {short} vs {a.ratio}x)")

for name, d in [("fg", fg), ("bg", bgm)]:
    out = d[["chr", "start", "end"]].copy()
    out["name"] = f"{name}_" + out["chr"] + "_" + out["start"].astype(str) + "_" + out["end"].astype(str)
    out["chr"] = out["chr"].str.removeprefix("chr")
    out.sort_values(["chr", "start"]).to_csv(os.path.join(a.out_dir, f"motif_{name}_{a.tag}.bed"),
                                             sep="\t", header=False, index=False)
summ = pd.DataFrame([{"set": n, "n": len(d), "len_mean": d["len"].mean(), "len_median": d["len"].median(),
                      "gc_mean": d["gc"].mean(), "gc_sd": d["gc"].std(), "cpg_per100_mean": d["cpg_per100"].mean()}
                     for n, d in [("foreground", fg), ("background_pool", bg), ("background", bgm)]])
summ.to_csv(os.path.join(a.out_dir, "background_match.tsv"), sep="\t", index=False, float_format="%.4g")
print(summ.to_string(index=False))
