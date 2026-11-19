#!/usr/bin/env bash
# One streaming pass over the CpG x TSS-window intersect (chr start end CpGID TxID distance): writes one
# TxID\tCpGID file per TxID (the candidate CpGs of that model), the TxID list and two small aggregate
# tables for the QC plots. Replaces filter_selected_TxIDs_second_round_NTM.R (fread of the whole table,
# ~600M rows / 33 GB for +-1 Mb windows -> killed for memory) and the per-gene grep of unpack_topK_NTM.
# Usage: split_topK_by_TxID.sh <intersect.bed> <split_dir> <ff_selected_TxIDs> <cpgs_per_txid.tsv> <distance_hist.tsv>
set -euo pipefail
export LC_ALL=C          # decimal point, not comma, in the histogram bins (the machine locale uses commas)
IN=$1; SPLIT=$2; TXIDS=$3; COUNTS=$4; HIST=$5
rm -rf "$SPLIT"; mkdir -p "$SPLIT" "$(dirname "$TXIDS")" "$(dirname "$COUNTS")" "$(dirname "$HIST")"
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true     # one output file per TxID stays open during the pass

mawk -F'\t' -v OFS='\t' -v dir="$SPLIT" -v cnt="$COUNTS" -v hist="$HIST" '
{
    tx = $5
    print tx, $4 > (dir "/" tx ".bed")
    n[tx]++
    d = $6 + 0; if (d < 1) d = 1
    b = int(log(d) / log(10) * 10)            # log10(distance) in bins of 0.1
    h[b]++
}
END {
    for (t in n) print t, n[t] > cnt
    for (b in h) print b / 10, h[b] > hist
}' "$IN"

cut -f1 "$COUNTS" | sort > "$TXIDS"
echo "$(wc -l < "$TXIDS") TxIDs, $(ls "$SPLIT" | wc -l) per-TxID files -> $SPLIT"
