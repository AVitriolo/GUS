import argparse

import pandas as pd

parser = argparse.ArgumentParser(description="Build PUREE input from the TMM+log2+ComBat-corrected "
                                             "transcript matrices, aggregated to gene-level HGNC symbols")
parser.add_argument("--combat_train", required=True)
parser.add_argument("--combat_test", required=True)
parser.add_argument("--conversion_table", required=True, help="transcript_id -> gene_name mapping TSV")
parser.add_argument("--sample_prefix", default=None,
                    help="if set, keep only sample columns starting with this prefix (e.g. 'N', 'T', 'DTB')")
parser.add_argument("--output", required=True)
args = parser.parse_args()

train = pd.read_csv(args.combat_train, sep="\t", index_col=0)
test = pd.read_csv(args.combat_test, sep="\t", index_col=0)
log_combat = pd.concat([train, test], axis=1)

if args.sample_prefix is not None:
    log_combat = log_combat.loc[:, log_combat.columns.str.startswith(args.sample_prefix)]
print(f"{log_combat.shape[1]} samples selected, {log_combat.shape[0]} transcripts")

# Undo log2(TMM-normalized + 1); ComBat's correction is not constrained to keep
# values positive on the log scale, so clip negative counts at 0 after inversion.
linear = (2 ** log_combat - 1).clip(lower=0)

conversion = pd.read_csv(args.conversion_table, sep=" ", usecols=["transcript_id", "gene_name"])
conversion = conversion.drop_duplicates(subset="transcript_id").set_index("transcript_id")

gene_name = conversion["gene_name"].reindex(linear.index)
matched = gene_name.notna()
print(f"{matched.sum()} / {len(linear)} transcripts mapped to a gene symbol")

gene_counts = linear.loc[matched].groupby(gene_name[matched]).sum()

gene_counts.T.to_csv(args.output, sep="\t", index_label="")
