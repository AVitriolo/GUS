#!/usr/bin/env bash
# Download the TCGA HumanMethylation450 and HiSeqV2 files of all TCGA cohorts from the UCSC Xena tcgaHub, on a
# machine WITH internet, count the samples having both, and write the cohort selection. Then copy OUT_DIR to the
# cluster's WORK_DIR (the prepare jobs unpack the .gz files and skip the download).
# Step 1 reads only the header line of every file (a few MB each) to count the samples with both data types;
# step 2 downloads the full files of the TOP_N cohorts with the most matched samples (and >= MIN_SAMPLES).
# Usage: ./download_local.sh <OUT_DIR> [TOP_N=15] [MIN_SAMPLES=100] [PARALLEL=4]
set -uo pipefail
OUT=${1:?usage: download_local.sh OUT_DIR [TOP_N] [MIN_SAMPLES] [PARALLEL]}; TOP=${2:-15}; MIN=${3:-100}; PAR=${4:-4}
COHORTS="ACC BLCA BRCA CESC CHOL COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM"
URL=https://tcga.xenahubs.net/download
mkdir -p "$OUT"
hdr() { curl -sL -m 120 "$URL/TCGA.$1.sampleMap/$2.gz" 2>/dev/null | head -c 8000000 | gunzip 2>/dev/null | head -1 | tr '\t' '\n' | tail -n +2 | sort -u; }
echo -e "cohort\tn_methylation\tn_expression\tn_both" > "$OUT/cohorts_selected.txt.tsv"
echo "[$(date +%T)] step 1: counting samples from the file headers"
for C in $COHORTS; do
  m=$(hdr $C HumanMethylation450); e=$(hdr $C HiSeqV2)
  nm=$(echo -n "$m" | grep -c .); ne=$(echo -n "$e" | grep -c .); nb=$(comm -12 <(echo "$m") <(echo "$e") | grep -c .)
  echo -e "$C\t$nm\t$ne\t$nb" >> "$OUT/cohorts_selected.txt.tsv"; echo "[$(date +%T)] $C methylation $nm | expression $ne | both $nb"
done
{ head -1 "$OUT/cohorts_selected.txt.tsv"; tail -n +2 "$OUT/cohorts_selected.txt.tsv" | sort -t$'\t' -k4,4nr; } > "$OUT/tmp.tsv" && mv "$OUT/tmp.tsv" "$OUT/cohorts_selected.txt.tsv"
awk -F'\t' -v m=$MIN 'NR>1 && $4>=m {print $1}' "$OUT/cohorts_selected.txt.tsv" | head -n "$TOP" > "$OUT/cohorts_selected.txt"
column -t "$OUT/cohorts_selected.txt.tsv"
echo "[$(date +%T)] selected (top $TOP with >= $MIN samples having both): $(tr '\n' ' ' < "$OUT/cohorts_selected.txt")"
get() {   # one cohort: download both files (resumable, retried), skip if already complete
  C=$1; D=$OUT/$C/TCGA.$C.sampleMap; mkdir -p "$D"
  for f in HumanMethylation450 HiSeqV2; do
    t=$D/$f.gz
    if [ -s "$t" ] && gzip -t "$t" 2>/dev/null; then echo "[$(date +%T)] $C $f: already downloaded"; continue; fi
    for try in 1 2 3; do
      curl -sL --retry 3 -o "$t" "$URL/TCGA.$C.sampleMap/$f.gz" && gzip -t "$t" 2>/dev/null && break
      echo "[$(date +%T)] $C $f: attempt $try failed, retrying"; rm -f "$t"; sleep 10
    done
    [ -s "$t" ] && gzip -t "$t" 2>/dev/null && echo "[$(date +%T)] $C $f: $(du -h "$t" | cut -f1)" || echo "[$(date +%T)] $C $f: FAILED"
  done
}
export -f get; export OUT URL
echo "[$(date +%T)] step 2: downloading the selected cohorts, $PAR at a time"
xargs -P "$PAR" -I{} bash -c 'get {}' < "$OUT/cohorts_selected.txt"
echo "Total size: $(du -sh "$OUT" | cut -f1). Copy $OUT/ to the cluster WORK_DIR, then: sbatch submit_selected.sbatch"
