"""CV (or train) R² vs. held-out test R² per transcript model.

--x_metric=cv_r2 (default) shows whether the cross-validation score used to
compare configurations (plot_tcga_prad_vs_wgbs_cv_r2.py) ranks models the same
way as the held-out test set. --x_metric=train_r2 shows the in-sample fit of
the final refit model against the test set, i.e. how much the models overfit.
With ~27k models a plain scatter saturates around 0, so points are drawn as a
log-scaled hexbin density; correlations are computed on all models, including
those outside the clipped axes.

Usage (2 Mb k-capped main run):
  python analyses/scripts/plot_cv_vs_test_r2.py \
    --performance_dir=results/performance_2mb_full_kcapped_cv3r5 \
    --prefix=hg38_v38_pc_50000_25_10_3 \
    --output_path=results/plots/wgbs_2mb_kcapped_cv_vs_test_r2.pdf
  python analyses/scripts/plot_cv_vs_test_r2.py --x_metric=train_r2 \
    --performance_dir=results/performance_2mb_full_kcapped_cv3r5 \
    --prefix=hg38_v38_pc_50000_25_10_3 \
    --output_path=results/plots/wgbs_2mb_kcapped_train_vs_test_r2.pdf
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import pearsonr, spearmanr

# One-hue light->dark ramp on the WGBS aqua (validated categorical slot 3),
# so density reads as magnitude and matches the WGBS boxes of Fig. 1B.
AQUA_RAMP = LinearSegmentedColormap.from_list("aqua", ["#d4f0e5", "#1baf7a", "#0a4a33"])
INK = "#0b0b0b"
MUTED = "#898781"
X_METRICS = {
    "cv_r2": ("CV", "CV R² (5× repeated 3-fold CV, training set)"),
    "train_r2": ("Train", "Train R² (final model, in-sample on the training set)"),
}


def read_metrics(path, keys):
    values = {}
    with open(path) as fh:
        for line in fh:
            key, _, value = line.strip().partition("\t")
            if key in keys:
                try:
                    values[key] = float(value)
                except ValueError:
                    values[key] = np.nan
    return values


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--performance_dir", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--x_metric", default="cv_r2", choices=list(X_METRICS))
    p.add_argument("--output_path", required=True)
    p.add_argument("--no_titles", action="store_true", help="omit the title (report version: it goes in the caption)")
    p.add_argument("--unit", default="transcript", help="model unit in labels (transcript / gene)")
    args = p.parse_args()
    xm = args.x_metric
    short, xlabel = X_METRICS[xm]

    files = sorted(glob.glob(os.path.join(args.performance_dir, f"{args.prefix}_*")))
    df = pd.DataFrame([read_metrics(f, (xm, "test_r2")) for f in files]).dropna()
    if df.empty:
        print(f"No {xm}/test_r2 pairs found in {args.performance_dir} - skipping plot.")
        return
    x, y = df[xm], df.test_r2

    rho = spearmanr(x, y).statistic
    r = pearsonr(x, y).statistic
    diff = x - y
    stats = f"Spearman ρ = {rho:.2f}\nPearson r = {r:.2f}\n"
    # the gap lines only matter for the overfitting (train) view
    if xm == "train_r2":
        stats += (f"median {short} − test = {diff.median():+.3f}\n"
                  f"{short} > test: {100 * (diff > 0).mean():.0f}% of models\n")
    stats += f"n = {len(df):,} models"

    # Both axes share one range so y = x is the diagonal; the long negative
    # tail of poorly-fit models is clipped at the 1st percentile.
    lo = min(x.quantile(0.01), y.quantile(0.01))
    hi = min(1.0, max(x.max(), y.max())) + 0.02
    inside = x.between(lo, hi) & y.between(lo, hi)
    n_outside = int((~inside).sum())

    os.makedirs(os.path.dirname(args.output_path), exist_ok=True)
    # Drawn at roughly report size with the same 8/9 pt text as the coherence
    # figure (plot_wgbs_vs_tcga_feature_coherence.py), so both print alike at 100%.
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8})
    # constrained layout: tight_layout + the equal aspect pushed the y label off
    # the canvas once the title was dropped
    fig, ax = plt.subplots(figsize=(4.6, 4.0), layout="constrained")
    hb = ax.hexbin(x[inside], y[inside], gridsize=60, extent=(lo, hi, lo, hi),
                   bins="log", mincnt=1, cmap=AQUA_RAMP, linewidths=0)
    cbar = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(f"{args.unit.capitalize()} models per bin (log scale)")
    cbar.outline.set_visible(False)

    ax.plot([lo, hi], [lo, hi], color=MUTED, linestyle="--", linewidth=1.2, zorder=3)
    ax.axhline(0, color=MUTED, linewidth=0.6, alpha=0.5, zorder=0)
    ax.axvline(0, color=MUTED, linewidth=0.6, alpha=0.5, zorder=0)
    # test R² above the x metric is rare in both views, so the upper-left
    # corner stays empty. The train view keys the dashed line as the legend
    # under the stats; the CV view leaves it to the caption.
    if xm == "train_r2":
        ax.lines[0].set_label("y = x")
        leg = ax.legend(frameon=False, loc="upper left", title=stats, title_fontsize=8)
        leg._legend_box.align = "left"
    else:
        ax.text(0.03, 0.97, stats, transform=ax.transAxes, ha="left", va="top", color=INK)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Test R² (held-out 30%)")
    if not args.no_titles:
        ax.set_title(f"{short} vs. held-out test R² per {args.unit} model\n"
                     f"axes clipped at 1st percentile; {n_outside} models outside, included in ρ")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.savefig(args.output_path, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output_path}")
    print(f"n = {len(df)}, Spearman rho = {rho:.4f}, Pearson r = {r:.4f}, "
          f"median {xm} = {x.median():.4f}, median test = {y.median():.4f}, "
          f"median {short} - test = {diff.median():.4f}, IQR [{diff.quantile(.25):.4f}, {diff.quantile(.75):.4f}], "
          f"{short} > test = {100 * (diff > 0).mean():.1f}%, outside axes = {n_outside}")


if __name__ == "__main__":
    main()
