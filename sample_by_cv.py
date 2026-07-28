import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", default="cv_r2_by_transcript.tsv", type=str)
parser.add_argument("--output", default="sampled_350_transcripts.txt", type=str)
parser.add_argument("--n_per_bin", default=70, type=int)
parser.add_argument("--seed", default=48047510, type=int)
parser.add_argument("--mode", choices=["fixed", "quantile"], default="fixed",
                     help="'fixed' = the four cv_r2 ranges below; 'quantile' = N equal-sized "
                          "groups by rank (guarantees n_per_bin per group)")
args = parser.parse_args()

df = pd.read_csv(args.input, sep=r"\s+", header=None, names=["tx_id", "cv_r2"])
df["cv_r2"] = pd.to_numeric(df["cv_r2"], errors="coerce")
df = df.dropna(subset=["cv_r2"])

if args.mode == "fixed":
    bins = [-0.4, -0.2, 0.0, 0.2, 0.4]
    labels = ["-0.4_to_-0.2", "-0.2_to_0", "0_to_0.2", "0.2_to_0.4"]
    df["bin"] = pd.cut(df["cv_r2"], bins=bins, labels=labels, right=False)
else:
    n_needed = 5 * args.n_per_bin
    if len(df) < n_needed:
        print(f"WARNING: only {len(df)} scored transcripts total, need {n_needed} "
              f"for {args.n_per_bin} per bin x 5 bins -- every bin will be short")
    df["rank_pct"] = df["cv_r2"].rank(method="first", pct=True)
    df["bin"] = pd.qcut(df["rank_pct"], q=5, labels=False, duplicates="drop")
    edges = df.groupby("bin", observed=True)["cv_r2"].agg(["min", "max"])
    labels = [f"q{int(b)} [{row['min']:.3f}-{row['max']:.3f}]" for b, row in edges.iterrows()]
    df["bin"] = df["bin"].map(dict(zip(edges.index, labels)))

rng = np.random.RandomState(args.seed)
sampled = []
for label in labels:
    pool = df[df["bin"] == label]
    print(f"{label}: {len(pool)} available")
    if len(pool) < args.n_per_bin:
        print(f"  WARNING: only {len(pool)} transcripts in bin '{label}', "
              f"taking all of them instead of {args.n_per_bin}")
        chosen = pool
    else:
        chosen = pool.sample(n=args.n_per_bin, random_state=rng)
    sampled.append(chosen)

result = pd.concat(sampled)
print(f"\nTotal sampled: {len(result)}")

with open(args.output, "w") as fh:
    for tx in result["tx_id"]:
        fh.write(f"{tx}\n")

result.to_csv(args.output.replace(".txt", "_with_r2.tsv"), sep="\t", index=False)
print(f"Wrote {args.output} and {args.output.replace('.txt', '_with_r2.tsv')}")
