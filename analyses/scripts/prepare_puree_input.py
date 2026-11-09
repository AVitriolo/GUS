import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Aggregate transcript-level kallisto counts to "
                                             "gene-level HGNC symbols, in PUREE's [samples, genes] orientation")
parser.add_argument("--counts", required=True, help="kallisto counts TSV: transcript_id rows x sample columns")
parser.add_argument("--conversion_table", required=True, help="transcript_id -> gene_name mapping TSV")
parser.add_argument("--sample_prefix", default=None,
                    help="if set, keep only sample columns starting with this prefix (e.g. 'N', 'T')")
parser.add_argument("--output", required=True)
args = parser.parse_args()

counts = pd.read_csv(args.counts, sep="\t", index_col="transcript_id")

if args.sample_prefix is not None:
    counts = counts.loc[:, counts.columns.str.startswith(args.sample_prefix)]
print(f"{counts.shape[1]} samples selected")

conversion = pd.read_csv(args.conversion_table, sep=" ", usecols=["transcript_id", "gene_name"])
conversion = conversion.drop_duplicates(subset="transcript_id").set_index("transcript_id")

gene_name = conversion["gene_name"].reindex(counts.index)
matched = gene_name.notna()
print(f"{matched.sum()} / {len(counts)} transcripts mapped to a gene symbol")

gene_counts = counts.loc[matched].groupby(gene_name[matched]).sum()

gene_counts.T.to_csv(args.output, sep="\t", index_label="")
