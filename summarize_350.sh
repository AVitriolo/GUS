#!/bin/bash
# Run from Prostate_GAS root (same directory as results/performance/)
set -euo pipefail

TX_LIST="diagnostics/sampled_350_transcripts.txt"
SRC_DIR="results/performance"
OUT="diagnostics/sampled_350_r2_summary.tsv"

echo -e "tx_id\ttrain_r2\tcv_r2\ttest_r2" > "$OUT"

n_ok=0
n_missing=0

while read -r TXID; do
    [ -z "$TXID" ] && continue

    matches=$(find "$SRC_DIR" -maxdepth 1 -name "*_${TXID}")
    n_match=$(echo "$matches" | grep -c . || true)

    if [ "$n_match" -eq 0 ]; then
        echo "MISSING: $TXID"
        echo -e "${TXID}\tNA\tNA\tNA" >> "$OUT"
        n_missing=$((n_missing + 1))
        continue
    fi
    if [ "$n_match" -gt 1 ]; then
        echo "AMBIGUOUS: $TXID matched multiple files, using first:"
        echo "$matches"
    fi

    FILE=$(echo "$matches" | head -n1)

    TRAIN_R2=$(awk -F'\t' '$1=="train_r2"{print $2}' "$FILE")
    CV_R2=$(awk -F'\t' '$1=="cv_r2"{print $2}' "$FILE")
    TEST_R2=$(awk -F'\t' '$1=="test_r2"{print $2}' "$FILE")

    TRAIN_R2=${TRAIN_R2:-NA}
    CV_R2=${CV_R2:-NA}
    TEST_R2=${TEST_R2:-NA}

    echo -e "${TXID}\t${TRAIN_R2}\t${CV_R2}\t${TEST_R2}" >> "$OUT"
    n_ok=$((n_ok + 1))
done < "$TX_LIST"

echo ""
echo "Done: $n_ok found, $n_missing missing. Wrote $OUT"
