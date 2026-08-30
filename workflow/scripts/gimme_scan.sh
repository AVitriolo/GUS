#!/usr/bin/env bash
# =============================================================================
# gimme_scan_metacpgs.sh
# TF motif scan on metaCpG regions.
#   0. (optional) restrict to a subset of cluster IDs (--subset file, one ID/line)
#   1. pad any region < MIN_LEN bp out to MIN_LEN, centered on its midpoint
#   2. strip the "chr" prefix (GRCh38.p14 FASTA uses Ensembl names: 1,2,..)
#   3. gimme scan against the genome
#
# Args (positional; use "" to skip an optional one):
#   $1  input_bed      region BED (chr start end name ...)
#   $2  genome         genomepy genome name (GRCh38.p14) or FASTA path
#   $3  output_tsv     gimme scan output
#   $4  chrom_sizes    genome .chrom.sizes (chr-prefixed; for the pad clamp)
#   $5  min_len        minimum region length (default 20)
#   $6  pfm            motif database (optional; gimme vertebrate v5 default)
#   $7  subset         OPTIONAL file of cluster IDs (BED col 4) to keep; "" = all
# =============================================================================
set -euo pipefail

input_bed="$1"
genome="$2"
output_tsv="$3"
chrom_sizes="$4"
min_len="${5:-20}"
pfm="${6:-}"
subset="${7:-}"

workdir="$(dirname "$output_tsv")"
mkdir -p "$workdir"
padded_bed="${output_tsv%.tsv}_padded.bed"
scan_input="$input_bed"

# ---- 0. optional subset by cluster ID (BED column 4) ----
if [ -n "$subset" ]; then
    if [ ! -s "$subset" ]; then
        echo "ERROR: subset file '$subset' is empty or missing" >&2
        exit 1
    fi
    subset_bed="${output_tsv%.tsv}_subset.bed"
    # keep BED rows whose name (col 4) is in the subset list
    awk 'NR==FNR { keep[$1]=1; next } ($4 in keep)' "$subset" "$input_bed" > "$subset_bed"
    echo "subset: kept $(wc -l < "$subset_bed") of $(wc -l < "$input_bed") regions"
    scan_input="$subset_bed"
fi

# ---- 1. pad (chr-prefixed, via chrom.sizes)  +  2. strip chr for Ensembl FASTA ----
awk -v OFS='\t' -v minlen="$min_len" '
    NR==FNR { clen[$1]=$2; next }
    {
        len = $3 - $2
        if (len < minlen) {
            mid  = int(($2 + $3) / 2)
            half = int(minlen / 2)
            s = mid - half
            e = mid + (minlen - half)
            if (s < 0) { e = e - s; s = 0 }
            cl = clen[$1]
            if (cl != "" && e > cl) { s = s - (e - cl); e = cl; if (s < 0) s = 0 }
            $2 = s; $3 = e
        }
        sub(/^chr/, "", $1)
        print
    }
' "$chrom_sizes" "$scan_input" > "$padded_bed"

echo "padded + renamed: $(wc -l < "$padded_bed") regions -> $padded_bed"

# ---- 3. gimme scan ----
if [ -n "$pfm" ]; then
    gimme scan "$padded_bed" -g "$genome" -p "$pfm" -b > "$output_tsv"
else
    gimme scan "$padded_bed" -g "$genome" -b > "$output_tsv"
fi

echo "gimme scan output: $output_tsv"
