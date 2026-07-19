#!/usr/bin/env python
# ============================================================
# split_samples.py
#
# Defines the ONE train/test partition used by BOTH modalities.
#
#   1. load raw RNA counts (N+T and M) and the WGBS colData
#   2. keep samples present in BOTH modalities  (paired samples only)
#   3. drop samples with RNA library size < min_libsize   (per-sample,
#      leakage-neutral, applied BEFORE the split)
#   4. split the survivors into train / test (deterministic, optionally
#      stratified by condition)
#   5. write train/test sample lists used by every downstream step
#
# Every parameter-learning step downstream (gene filter, TMM reference,
# CpG coverage filter, imputation medians, batch correction) is fit on
# the TRAIN list and applied to the TEST list.
# ============================================================
import argparse
import os

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

p = argparse.ArgumentParser()
p.add_argument("--input_expr_NT", required=True,
               help="Raw RNA counts TSV (N+T): transcript_id column + sample columns")
p.add_argument("--input_expr_M", required=True,
               help="Raw RNA counts TSV (M)")
p.add_argument("--input_wgbs_samples", required=True,
               help="Text file with one WGBS sample ID per line "
                    "(exported from integrated_wgbs_colData.rds)")
p.add_argument("--min_libsize", type=float, default=5e6)
p.add_argument("--test_size", type=float, default=0.3)
p.add_argument("--seed", type=int, default=48047510)
p.add_argument("--stratify", action="store_true",
               help="Stratify the split by condition (recommended)")
p.add_argument("--outdir", default="resources/split")
p.add_argument("--out_plot", default=None)
args = p.parse_args()

os.makedirs(args.outdir, exist_ok=True)


def load_clean(path):
    df = pd.read_csv(path, sep="\t")
    df["transcript_id"] = df["transcript_id"].str.replace(r"\.[0-9]*", "",
                                                          regex=True)
    return df.set_index("transcript_id")


def cond(s):
    if s.startswith("N"):
        return "Normal"
    if s.startswith("T"):
        return "Tumour"
    if s.startswith("DTB"):
        return "metastasis"
    return "unknown"


CONDS = ["Normal", "Tumour", "metastasis"]

# ---- 1. load ------------------------------------------------------------
NT = load_clean(args.input_expr_NT)
M = load_clean(args.input_expr_M)
rna = NT.join(M, how="outer")
rna_samples = list(rna.columns)

with open(args.input_wgbs_samples) as f:
    wgbs_samples = [l.strip() for l in f if l.strip()]

print(f"RNA samples : {len(rna_samples)}")
print(f"WGBS samples: {len(wgbs_samples)}")

# ---- 2. PAIRED samples only (present in both modalities) ----------------
paired = sorted(set(rna_samples) & set(wgbs_samples))
print(f"Paired (RNA & WGBS): {len(paired)}")
print(f"  RNA-only dropped : {len(set(rna_samples) - set(wgbs_samples))}")
print(f"  WGBS-only dropped: {len(set(wgbs_samples) - set(rna_samples))}")

if not paired:
    raise ValueError("No samples present in both RNA and WGBS.")

# ---- 3. library-size filter (per-sample; BEFORE the split) --------------
libsize = rna[paired].sum(axis=0)
keep = libsize[libsize >= args.min_libsize].index.tolist()
dropped = sorted(set(paired) - set(keep))
print(f"\nLibrary-size filter (>= {args.min_libsize:.0e}): "
      f"kept {len(keep)}, dropped {len(dropped)}")
if dropped:
    print("  dropped:", ", ".join(dropped))

conditions = [cond(s) for s in keep]
if "unknown" in conditions:
    unknown = [s for s, c in zip(keep, conditions) if c == "unknown"]
    raise ValueError(f"Could not infer condition for: {unknown[:10]}")

# ---- 4. deterministic split --------------------------------------------
train_samples, test_samples = train_test_split(
    keep,
    test_size=args.test_size,
    random_state=args.seed,
    stratify=conditions if args.stratify else None,
)
train_samples = sorted(train_samples)
test_samples = sorted(test_samples)

# ---- 5. write the lists -------------------------------------------------
for name, lst in [("train", train_samples),
                  ("test", test_samples),
                  ("kept", keep),
                  ("paired", paired)]:
    with open(os.path.join(args.outdir, f"{name}_samples.txt"), "w") as f:
        f.write("\n".join(lst) + "\n")


def bal(lst):
    c = [cond(s) for s in lst]
    return {k: c.count(k) for k in CONDS}


tr, te = bal(train_samples), bal(test_samples)
print(f"\nSplit (seed={args.seed}, test_size={args.test_size}, "
      f"stratify={args.stratify}):")
print(f"  Train: {len(train_samples):4d}  {tr}")
print(f"  Test : {len(test_samples):4d}  {te}")

pd.DataFrame({
    "split": ["train", "test"],
    "n_total": [len(train_samples), len(test_samples)],
    **{c: [tr[c], te[c]] for c in CONDS},
}).to_csv(os.path.join(args.outdir, "split_summary.tsv"),
          sep="\t", index=False)

# ---- optional QC plot ---------------------------------------------------
if args.out_plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(args.out_plot) or ".", exist_ok=True)

    x = np.arange(2)
    width = 0.25
    cols = {"Normal": "#2166ac", "Tumour": "#d6604d", "metastasis": "#1a9641"}

    fig, ax = plt.subplots(figsize=(8, 5))
    for i, c in enumerate(CONDS):
        bars = ax.bar(x + (i - 1) * width, [tr[c], te[c]], width,
                      label=c, color=cols[c])
        ax.bar_label(bars, fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([f"Train (n={len(train_samples)})",
                        f"Test (n={len(test_samples)})"])
    ax.set_ylabel("Number of samples")
    ax.set_title(f"Paired RNA+WGBS samples: train/test composition\n"
                 f"(seed={args.seed}, test_size={args.test_size}, "
                 f"stratify={args.stratify})")
    ax.legend(title="Condition")
    plt.tight_layout()
    plt.savefig(args.out_plot)
    print(f"  plot  -> {args.out_plot}")

print(f"\nSaved to {args.outdir}/: "
      f"train_samples.txt, test_samples.txt, kept_samples.txt, "
      f"paired_samples.txt, split_summary.tsv")