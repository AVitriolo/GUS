"""Stability selection with the models' own selector: shared core for any cohort (TCGA 450K, GBM WGBS, ...).

Same procedure as the prostate WGBS analysis (Prostate_GAS analyses/scripts/stability_selection_kbest.py and
stability_regions_corr.py), cohort-independent:
  * half-samples: n_pairs complementary pairs of the TRAIN samples, each half stratified by a sample group
    (tumour / normal, subtype, ...) - Shah & Samworth complementary-pairs stability selection
  * selector on each half: SelectKBest(f_regression) with the gene's tuned k (= number of its saved selected
    features), i.e. the selection step of model_selection_split_kbest_repeated_kcapped.py
  * co-methylated regions: candidates sorted by position; a new region starts when the next feature is more than
    max_gap bp away or its Pearson r with the previous one (train values) is < corr; expression is not used
Writes, per shard, the files read by stability_selection_summary.R and stability_selection_regions.R:
  subsample_selections_shard<S>.tsv  gene, half index, comma-separated 0-based candidate indices
  candidates_shard<S>.tsv            gene, n candidates, k, index, feature ID, saved_selected (features selected in
                                     >= 1 half or saved by the model only)
  region_map_corr<r>_shard<S>.tsv    gene, n regions, feature ID, region, region start, region end, n features, chr
"""
import zlib

import numpy as np
from sklearn.feature_selection import SelectKBest, f_regression

SEED = 48047510


def complementary_halves(groups, n_pairs, seed):
    """n_pairs x 2 index arrays: each pair splits every group in two halves (sizes differ by <= 1)."""
    groups = np.asarray(groups)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_pairs):
        h = np.zeros(len(groups), dtype=bool)
        for g in np.unique(groups):
            idx = rng.permutation(np.flatnonzero(groups == g))
            h[idx[: len(idx) // 2]] = True
        out += [np.flatnonzero(h), np.flatnonzero(~h)]
    return out


def select_k(X, y, k):
    return SelectKBest(f_regression, k=min(k, X.shape[1])).fit(X, y).get_support(indices=True)


def comethylated_regions(X, chrom, pos, max_gap=300, corr=0.5):
    """Region id (1-based) per column of X; columns must be sorted by (chrom, pos)."""
    Xc = X - X.mean(0)
    sd = np.sqrt((Xc ** 2).sum(0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (Xc[:, :-1] * Xc[:, 1:]).sum(0) / (sd[:-1] * sd[1:])
    chrom = np.asarray(chrom)
    brk = (chrom[1:] != chrom[:-1]) | (np.diff(pos) > max_gap) | ~(r >= corr)      # NaN r (constant) -> break
    return np.cumsum(np.concatenate([[True], brk]))


class ShardWriter:
    def __init__(self, out_dir, shard, corr=0.5):
        self.f_sel = open(f"{out_dir}/subsample_selections_shard{shard}.tsv", "w")
        self.f_cand = open(f"{out_dir}/candidates_shard{shard}.tsv", "w")
        self.f_map = open(f"{out_dir}/region_map_corr{corr}_shard{shard}.tsv", "w")

    def gene(self, gene, ids, chrom, pos, X, y, saved, groups, n_pairs, max_gap=300, corr=0.5):
        """ids/chrom/pos: candidate features sorted by position; X: train samples x candidates; saved: set of IDs."""
        k = len(saved)
        if k == 0 or X.shape[1] == 0:
            return False
        halves = complementary_halves(groups, n_pairs, [SEED, zlib.crc32(gene.encode())])
        ever = set()
        for b, h in enumerate(halves):
            s = select_k(X[h], y[h], k)
            ever.update(s.tolist())
            self.f_sel.write(f"{gene}\t{b}\t{','.join(map(str, s))}\n")
        keep = [n for n, c in enumerate(ids) if n in ever or c in saved]
        for n in keep:
            self.f_cand.write(f"{gene}\t{len(ids)}\t{k}\t{n}\t{ids[n]}\t{int(ids[n] in saved)}\n")
        region = comethylated_regions(X, chrom, pos, max_gap, corr)
        starts = np.flatnonzero(np.concatenate([[True], region[1:] != region[:-1]]))
        ends = np.append(starts[1:], len(region)) - 1
        for n in keep:
            rg = region[n] - 1
            self.f_map.write(f"{gene}\t{region[-1]}\t{ids[n]}\t{region[n]}\t{pos[starts[rg]]}\t{pos[ends[rg]]}\t"
                             f"{ends[rg] - starts[rg] + 1}\t{chrom[n]}\n")
        for f in (self.f_sel, self.f_cand, self.f_map):
            f.flush()
        return True

    def close(self):
        for f in (self.f_sel, self.f_cand, self.f_map):
            f.close()
