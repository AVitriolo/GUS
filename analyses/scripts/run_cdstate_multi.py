# CDState multi-run (k range x restarts), following notebooks/CDState_end_to_end.ipynb
# (KMeans-on-cosine initialisation, Step 1 then purity-constrained Step 2).
# Each (k, iter) is independent and checkpointed: existing outputs are skipped.
import os, sys, argparse, copy, random, time
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[v] = "1"
os.environ["XLA_FLAGS"] = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"
sys.path.insert(0, "/home/Prostate_GAS/analyses/tools/CDState")
import numpy as np, pandas as pd
from joblib import Parallel, delayed

p = argparse.ArgumentParser()
p.add_argument("--bulk"); p.add_argument("--purity"); p.add_argument("--outdir"); p.add_argument("--prefix", default="M")
p.add_argument("--kmin", type=int, default=3); p.add_argument("--kmax", type=int, default=8)
p.add_argument("--niter", type=int, default=10); p.add_argument("--ncore", type=int, default=6)
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)

ann = pd.read_csv("/home/Prostate_GAS/analyses/tools/CDState/data/gene_order.csv")
df = pd.read_csv(a.bulk, index_col=0)
df = df.loc[df.index.isin(ann.gene_name)]
purity = pd.read_csv(a.purity, index_col=0)["purity"].loc[df.columns]
purity.index = df.columns
print(f"input: {df.shape[0]} annotated genes x {df.shape[1]} samples", flush=True)

def one(k, i):
    import CDState_base as cd
    from sklearn.metrics.pairwise import cosine_similarity
    from sklearn.cluster import KMeans
    stem = os.path.join(a.outdir, f"{a.prefix}_k{k}_iter{i}")
    if all(os.path.exists(f"{stem}_{s}.csv") for s in ("proportions", "sources", "mal")): return (k, i, "skip", 0)
    t0 = time.time(); random.seed(1000 * k + i); np.random.seed(1000 * k + i)
    cn = cd.CDState(df, num_bases=k, global_round=False); cn.prepare_data()
    lab = KMeans(n_clusters=k, n_init="auto", init="random", random_state=1000 * k + i).fit(cosine_similarity(cn.data.T)).labels_
    idx = [random.sample(np.where(lab == y)[0].tolist(), 1)[0] for y in range(k)]
    cn.W = copy.copy(cn.data[:, idx]) + 1e-10; cn.factorize()
    cg = cd.CDState(df, purity, global_round=True, num_bases=k, gene_list=cn.gene_list); cg.prepare_data()
    cg.W = copy.copy(cn.W) + 1e-10; cg.H = copy.copy(cn.H); cg.factorize()
    pd.DataFrame(cg.W, index=cg.gene_list).to_csv(f"{stem}_sources.csv.tmp")
    np.savetxt(f"{stem}_proportions.csv.tmp", cg.H, delimiter=",", header=",".join(cg.samples))
    np.savetxt(f"{stem}_mal.csv.tmp", np.atleast_1d(cg.mal), delimiter=",")
    for s in ("sources", "proportions", "mal"): os.replace(f"{stem}_{s}.csv.tmp", f"{stem}_{s}.csv")
    return (k, i, "done", time.time() - t0)

jobs = [(k, i) for k in range(a.kmin, a.kmax + 1) for i in range(a.niter)]
for k, i, st, dt in Parallel(n_jobs=a.ncore, verbose=0, return_as="generator_unordered")(delayed(one)(k, i) for k, i in jobs):
    print(f"k={k} iter={i} {st} {dt/60:.1f} min", flush=True)
print("all runs finished", flush=True)
