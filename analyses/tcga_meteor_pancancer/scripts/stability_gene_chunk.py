"""Stability selection for one chunk of genes of one cohort (one task of slurm/stability.sbatch).

The per-gene train tables are deleted after each model (KEEP_GENE_DATASETS=0), so they are rebuilt with
03_make_gene_datasets.py, the same construction and seed as the models (identical train samples and candidates).
Then, per gene with a model (results_<scale>/sel_feats/<gene>_SVR.tsv), stability_core: STAB_N_PAIRS complementary
pairs of half-samples of the train samples, stratified by sample type from the TCGA barcode (01-09 tumour, 10-19
normal, 20+ control), SelectKBest(f_regression) with the model's k on each half, co-methylated regions of the
candidate probes (hg19 positions; <= STAB_MAX_GAP bp, adjacent r >= STAB_CORR).
Writes <work_dir>/stability_<scale>/{subsample_selections,candidates,region_map_corr<r>}_shard<chunk>.tsv and a
done_<chunk> marker (a finished chunk is skipped when the task is rerun).
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stability_core import ShardWriter

p = argparse.ArgumentParser()
for k in ("cohort", "work_dir", "pipe_dir", "chunk"):
    p.add_argument(f"--{k}", required=True)
p.add_argument("--scale", default="beta")
p.add_argument("--window_bp", type=int, default=1_000_000)
p.add_argument("--test_size", type=float, default=0.3)
p.add_argument("--n_pairs", type=int, default=100)
p.add_argument("--corr", default="0.5")
p.add_argument("--max_gap", type=int, default=300)
a = p.parse_args()

C, R = a.work_dir, f"{a.work_dir}/results_{a.scale}"
OUT = f"{a.work_dir}/stability_{a.scale}"
shard = os.path.basename(a.chunk).replace(".txt", "")
os.makedirs(OUT, exist_ok=True)
if os.path.exists(f"{OUT}/done_{shard}"):
    print(f"{shard}: already done"); sys.exit(0)

genes = [g.strip() for g in open(a.chunk) if g.strip()]
genes = [g for g in genes if os.path.exists(f"{R}/sel_feats/{g}_SVR.tsv")]
print(f"[{a.cohort} {shard}] {len(genes)} genes with a model", flush=True)
tmp = f"{OUT}/tmp_{shard}"
shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
with open(f"{tmp}/genes.txt", "w") as fh:
    fh.write("\n".join(genes) + "\n")
if genes:
    subprocess.run([sys.executable, f"{a.pipe_dir}/scripts/03_make_gene_datasets.py", f"--cohort={a.cohort}",
                    f"--data_dir={C}", f"--gene_assoc={a.pipe_dir}/resources/HM450.hg19.gene_association.tsv.gz",
                    f"--genes_file={tmp}/genes.txt", f"--out_dir={tmp}", f"--window_bp={a.window_bp}",
                    f"--test_size={a.test_size}", f"--scale={a.scale}"], check=True, stdout=subprocess.DEVNULL)

rng = pd.read_csv(f"{C}/{a.cohort}_data_meth_range.csv", usecols=["seqnames", "start", "probe.id"])
where = dict(zip(rng["probe.id"], zip(rng["seqnames"], rng["start"].astype(int))))


def sample_group(s):
    """TCGA barcode sample type: TCGA-XX-XXXX-01 -> tumour; 10-19 normal; 20+ control; otherwise one group."""
    try:
        code = int(s.split("-")[3][:2])
    except (IndexError, ValueError):
        return "all"
    return "tumour" if code < 10 else "normal" if code < 20 else "control"


w = ShardWriter(OUT, shard, a.corr)
n = 0
for g in genes:
    f = f"{tmp}/{g}_train.tsv"
    if not os.path.exists(f):
        continue
    tr = pd.read_csv(f, sep="\t", index_col=0)
    y = tr.pop(g).to_numpy(dtype=np.float64)
    probes = [c for c in tr.columns if c in where]
    probes.sort(key=lambda c: where[c])
    chrom = [where[c][0] for c in probes]
    pos = np.array([where[c][1] for c in probes])
    saved = {l.strip() for l in open(f"{R}/sel_feats/{g}_SVR.tsv") if l.strip()}
    groups = [sample_group(s) for s in tr.index]
    if w.gene(g, probes, chrom, pos, tr[probes].to_numpy(dtype=np.float64), y, saved, groups, a.n_pairs,
              a.max_gap, float(a.corr)):
        n += 1
w.close()
shutil.rmtree(tmp, ignore_errors=True)
open(f"{OUT}/done_{shard}", "w").write(f"{n}\n")
print(f"[{a.cohort} {shard}] done: {n} genes", flush=True)
