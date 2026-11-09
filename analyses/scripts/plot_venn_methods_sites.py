"""Bone vs lymph node MethylDriver genes (q < Q) in bulk / CAMDAC / PureBeta, all on non-ComBat betas with identical CpGs,
samples and regions (run_methyldriver_purebeta_camdac_nocombat.sh). Gene level: multi-gene regions split; a gene's
direction is the sign of its strongest significant region (Δ = LN − bone; > 0 higher in lymph node).
Writes <out_prefix>_membership.tsv (gene, region, in_bulk/in_camdac/in_purebeta, direction per method) and the Venn,
with the genes of small regions written in the figure.
Usage: python plot_venn_methods_sites.py --res_dir DIR --out_prefix PREFIX [--q 0.05] [--max_list 12]
"""
import argparse
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from matplotlib.colors import to_rgba

p = argparse.ArgumentParser()
p.add_argument("--res_dir", required=True); p.add_argument("--out_prefix", required=True)
p.add_argument("--q", type=float, default=.05); p.add_argument("--max_list", type=int, default=12)
p.add_argument("--tag", default="2mb_v2_beta_genestd_stable_corr0.5_w20")
p.add_argument("--comparison", default="D_Bone_vs_D_Lymph", help="D_Bone_vs_D_Lymph or N_vs_T")
a = p.parse_args()
M = {"bulk": "bulk", "camdac": "CAMDAC", "purebeta": "PureBeta"}
dirs = {}
for k in M:
    d = pd.read_csv(f"{a.res_dir}/{k}/full_{a.comparison}_{a.tag}.tsv", sep="\t")
    d = d[d.q_value < a.q].assign(gene=lambda x: x.gene.astype(str).str.split(",")).explode("gene")
    d = d[d.gene.ne("nan")]
    d = d.iloc[np.argsort(-d.home_region_methyl_change.abs().to_numpy(), kind="stable")].drop_duplicates("gene")
    dirs[k] = d.set_index("gene").home_region_methyl_change.map(lambda v: (("higher in lymph node" if v > 0 else "higher in bone") if a.comparison != "N_vs_T" else ("higher in tumour" if v > 0 else "lower in tumour")))
genes = sorted(set().union(*[set(v.index) for v in dirs.values()]))
T = pd.DataFrame({"gene": genes})
for k in M: T[f"in_{k}"] = T.gene.isin(dirs[k].index); T[f"dir_{k}"] = T.gene.map(dirs[k])
T["region"] = T.apply(lambda r: "+".join(M[k] for k in M if r[f"in_{k}"]) + (" only" if sum(r[f"in_{k}"] for k in M) == 1 else ""), axis=1)
T.to_csv(f"{a.out_prefix}_membership.tsv", sep="\t", index=False)
print(T.region.value_counts().to_string())
for reg, x in T.groupby("region"):
    dcol = [f"dir_{k}" for k in M if x[f"in_{k}"].all()]
    agree = (x[dcol].nunique(axis=1) == 1).mean() if len(dcol) > 1 else np.nan
    print(f"{reg}: {len(x)} genes | " + " | ".join(f"{c[4:]}: {dict(x[c].value_counts())}" for c in dcol) + (f" | same direction {100*agree:.0f}%" if len(dcol) > 1 else ""))

# ---- Venn: thin dark outlines, clear tints for each method's own area, cream pairwise overlaps, amber core ----
C = {"bulk": (-.8, .5), "camdac": (.8, .5), "purebeta": (0, -.85)}; R = 1.5
COL = {"bulk": "#3f8fd2", "camdac": "#e8775a", "purebeta": "#2fae8f"}       # clear blue / coral / teal
xx, yy = np.meshgrid(np.linspace(-2.6, 2.6, 1300), np.linspace(-2.95, 2.6, 1400))
ins = {k: (xx - x) ** 2 + (yy - y) ** 2 <= R ** 2 for k, (x, y) in C.items()}
img = np.ones(xx.shape + (4,)); img[..., 3] = 0
for k in M:
    only = ins[k] & ~np.logical_or.reduce([ins[j] for j in M if j != k]); img[only] = to_rgba(COL[k], .26)
for i, j in (("bulk", "camdac"), ("bulk", "purebeta"), ("camdac", "purebeta")):
    o = [x for x in M if x not in (i, j)][0]; img[ins[i] & ins[j] & ~ins[o]] = to_rgba("#f6e7c8", .9)
img[ins["bulk"] & ins["camdac"] & ins["purebeta"]] = to_rgba("#f2b84b", .85)
plt.rcParams.update({"font.size": 12, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(7.2, 7.2)); ax.imshow(img, extent=(-2.6, 2.6, -2.95, 2.6), origin="lower", zorder=0, interpolation="antialiased")
for k, (x, y) in C.items(): ax.add_patch(plt.Circle((x, y), R, fill=False, lw=1.0, color="#3a3a37"))
pos = {"bulk only": (-1.5, 1.0), "CAMDAC only": (1.5, 1.0), "PureBeta only": (0, -1.7), "bulk+CAMDAC": (0, 1.25),
       "bulk+PureBeta": (-.95, -.45), "CAMDAC+PureBeta": (.95, -.45), "bulk+CAMDAC+PureBeta": (0, .1)}
for reg, (x, y) in pos.items():
    n = int((T.region == reg).sum()); core = reg == "bulk+CAMDAC+PureBeta"
    ax.text(x, y, f"{n}", ha="center", va="center", fontsize=19 if core else 17, color="#1f1f1d", fontweight="bold" if core else "normal")
for k, (tx, ty) in zip(M, ((-1.75, 2.2), (1.75, 2.2), (0, -2.5))):
    ax.text(tx, ty, M[k], ha="center", va="bottom" if ty > 0 else "top", fontsize=15, color=COL[k], fontweight="bold")
    ax.text(tx, ty - (.02 if ty > 0 else .2), f"{T[f'in_{k}'].sum()} genes", ha="center", va="top", fontsize=11, color="#52514e")
ax.set_xlim(-2.6, 2.6); ax.set_ylim(-3.05, 2.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title(f"{'Normal vs tumour' if a.comparison == 'N_vs_T' else 'Bone vs lymph node'} · MethylDriver genes, q < {a.q:g}", fontsize=13, color="#1f1f1d", pad=4)
fig.tight_layout()
for e in ("png", "pdf"): fig.savefig(f"{a.out_prefix}.{e}", dpi=220, bbox_inches="tight")
