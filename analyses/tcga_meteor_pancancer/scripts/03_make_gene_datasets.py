"""Per-gene train/test tables for a chunk of genes (generalised make_datasets_prad_positional_cli.py).

Same construction as the TCGA PRAD analysis:
  * TSS anchor = the 450K probe with the smallest |distToTSS| for the gene in the Zhou-lab hg19 gene-association table
  * features = ALL probes of the cohort's methylation matrix within +-WINDOW_BP of that anchor (position based)
  * target = the gene's expression (Xena HiSeqV2, gene level)
  * samples = all samples with both data types, split 70/30 with the fixed seed 48047510
Differences: features are BETA values by default (--scale=M for the original M-values), and the methylation matrix is read once per chunk instead of grepping the CSV once per gene (same result, faster).

Usage: python 03_make_gene_datasets.py --cohort=BRCA --data_dir=<WORK_DIR>/BRCA --gene_assoc=<resources>/HM450.hg19.gene_association.tsv.gz
           --genes_file=chunk_000.txt --out_dir=<WORK_DIR>/BRCA/datasets --window_bp=1000000 --test_size=0.3
"""
import argparse, os
import pandas as pd
from sklearn.model_selection import train_test_split

p = argparse.ArgumentParser()
for k in ("cohort", "data_dir", "gene_assoc", "genes_file", "out_dir"): p.add_argument(f"--{k}", required=True)
p.add_argument("--window_bp", type=int, default=1_000_000); p.add_argument("--test_size", type=float, default=0.3)
p.add_argument("--scale", choices=["beta", "M"], default="beta", help="methylation values used as features")
a = p.parse_args()
SEED = 48047510
genes = [l.strip() for l in open(a.genes_file) if l.strip()]
os.makedirs(a.out_dir, exist_ok=True)

ga = pd.read_csv(a.gene_assoc, sep="\t", dtype=str)
cols = ["geneNames", "transcriptTypes", "transcriptIDs", "distToTSS"]
for c in cols: ga[c] = ga[c].str.split(";")
ex = ga.explode(cols).reset_index(drop=True).rename(columns={"geneNames": "geneName"})
ex["geneName"] = ex["geneName"].replace({"": pd.NA, "NA": pd.NA}); ex["distToTSS"] = pd.to_numeric(ex["distToTSS"], errors="coerce")
ex = ex.dropna(subset=["geneName", "distToTSS"]); ex = ex[ex.geneName.isin(genes)]

rng = pd.read_csv(f"{a.data_dir}/{a.cohort}_data_meth_range.csv"); rng.index = list(rng["Unnamed: 0"]); rng = rng.drop(columns="Unnamed: 0")
exp = pd.read_csv(f"{a.data_dir}/need_gene_exp_{a.cohort}.csv", index_col=0)

# probes needed by this chunk, then read only those rows of the methylation matrix
need, anchors = set(), {}
for g in genes:
    r = ex[ex.geneName == g]
    if len(r) == 0: print(f"[{g}] SKIPPED: gene not in the gene-association table (no TSS anchor)"); continue
    ar = r.loc[r.distToTSS.abs().idxmin()]; ch, pos = ar.CpG_chrm, int(ar.CpG_beg)
    sc = rng[rng.seqnames == ch]; inw = sc[(sc.iloc[:, 1].astype(int) >= pos - a.window_bp) & (sc.iloc[:, 1].astype(int) <= pos + a.window_bp)]
    probes = inw.iloc[:, 5].tolist(); anchors[g] = (ch, pos, probes); need.update(probes)
mfile = f"{a.data_dir}/{a.cohort}_data_meth_beta.csv" if a.scale == "beta" else f"{a.data_dir}/{a.cohort}_data_meth.csv"
it = pd.read_csv(mfile, index_col=0, chunksize=20000)
meth = pd.concat([c[c.index.isin(need)] for c in it]).astype("float32")
print(f"chunk: {len(genes)} genes, {len(anchors)} with an anchor, {meth.shape[0]:,} probes read")

for g, (ch, pos, probes) in anchors.items():
    m = meth.loc[meth.index.intersection(probes)]
    if m.shape[0] == 0: print(f"[{g}] SKIPPED: no probes within +-{a.window_bp} bp"); continue
    if g not in exp.index: print(f"[{g}] SKIPPED: not in the expression matrix"); continue
    y = exp.loc[g]
    if isinstance(y, pd.DataFrame): y = y.iloc[0]          # duplicated gene names after synonym search: keep the first row
    X = m.T; common = sorted(set(X.index) & set(y.index))
    full = pd.concat([X.loc[common], y.loc[common].rename(g)], axis=1)
    tr, te = train_test_split(full, test_size=a.test_size, random_state=SEED)
    tr.to_csv(f"{a.out_dir}/{g}_train.tsv", sep="\t"); te.to_csv(f"{a.out_dir}/{g}_test.tsv", sep="\t")
    print(f"[{g}] anchor {ch}:{pos} | {m.shape[0]} probes | train {tr.shape} test {te.shape}")
