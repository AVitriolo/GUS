"""Compare cv_r2 distributions: HM450-probe-restricted run vs. the original
full-CpG pipeline, for whichever TxIDs have completed so far in the HM450 run.

Overwrites a single fixed output file each call (no timestamping) so it always
reflects the latest progress - intended to be re-run periodically (e.g. every
100 completed TxIDs) by run_hm450_validation.sh, not accumulated as a series.

Usage:
  python analyses/scripts/plot_hm450_vs_original_cv_r2.py \
    --performance_dir_hm450=results/performance_hm450 \
    --performance_dir_original=results/performance \
    --prefix=hg38_v38_pc_50000_25_10 \
    --cv=3 \
    --output_path=results/plots/hm450_vs_original_cv_r2.png
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Okabe-Ito colorblind-safe pair.
COLOR_HM450 = "#E69F00"
COLOR_ORIGINAL = "#0072B2"


def read_cv_r2(path):
    with open(path) as fh:
        for line in fh:
            key, _, value = line.strip().partition("\t")
            if key == "cv_r2":
                return float(value)
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--performance_dir_hm450", required=True)
    p.add_argument("--performance_dir_original", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--cv", required=True)
    p.add_argument("--output_path", required=True)
    args = p.parse_args()

    pattern = os.path.join(args.performance_dir_hm450, f"{args.prefix}_{args.cv}_*")
    hm450_files = sorted(glob.glob(pattern))

    rows = []
    n_missing_original = 0
    for f in hm450_files:
        base = os.path.basename(f)
        tx_id = base[len(f"{args.prefix}_{args.cv}_"):]
        r2_hm450 = read_cv_r2(f)
        if r2_hm450 is None:
            continue

        original_path = os.path.join(args.performance_dir_original, base)
        if not os.path.exists(original_path):
            n_missing_original += 1
            continue
        r2_original = read_cv_r2(original_path)
        if r2_original is None:
            continue

        rows.append((tx_id, "HM450 probes", r2_hm450))
        rows.append((tx_id, "Original (full CpGs)", r2_original))

    if not rows:
        print("No matched TxIDs yet (need >=1 completed HM450 result with a "
              "matching original results/performance file) - skipping plot.")
        return

    import pandas as pd
    df = pd.DataFrame(rows, columns=["TxID", "pipeline", "cv_r2"])
    n_txids = df["TxID"].nunique()

    os.makedirs(os.path.dirname(args.output_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.boxplot(
        data=df, x="pipeline", y="cv_r2",
        order=["Original (full CpGs)", "HM450 probes"],
        palette={"Original (full CpGs)": COLOR_ORIGINAL, "HM450 probes": COLOR_HM450},
        width=0.5, fliersize=0, ax=ax,
    )
    sns.stripplot(
        data=df, x="pipeline", y="cv_r2",
        order=["Original (full CpGs)", "HM450 probes"],
        color="black", alpha=0.35, size=3, jitter=0.15, ax=ax,
    )
    ax.axhline(0, color="grey", linewidth=1, linestyle="--", zorder=0)
    ax.set_xlabel("")
    ax.set_ylabel("cv_r2")
    ax.set_title(f"cv_r2: HM450-restricted vs. original pipeline\n(n = {n_txids} TxIDs completed so far)")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(args.output_path, dpi=150)
    plt.close(fig)

    print(f"Written {args.output_path} (n={n_txids} matched TxIDs, "
          f"{n_missing_original} HM450 results had no matching original file)")


if __name__ == "__main__":
    main()
