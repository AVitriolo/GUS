#!/usr/bin/env bash
# =============================================================================
# bigwig_stream_average.sh
# Averages each ENCODE bigWig over a BED without keeping the bigWigs locally
# (~600 GB; they are only on the HPC): per file, download from the public
# ENCODE S3 bucket -> bigWigAverageOverBed -> keep name + mean (col 6, covered
# bases only, as make_feature_matrix_from_bigwigs_and_bed.R) -> delete.
# Resumable: accessions with an output are skipped; failed ones are retried on
# the next run.
# A local path instead of a URL (e.g. the HPC folder mounted with sshfs) is
# read in place: no download, and the file is never deleted.
#
# Args: $1 urls_tsv   accession <TAB> url or local path [<TAB> size]
#       $2 bed        regions, unique names in col 4
#       $3 out_dir    <accession>.tsv.gz (name, mean)
#       $4 tmp_dir    downloads (one bigWig per parallel job at a time)
#       $5 n_parallel
# bigWigAverageOverBed and curl must be on PATH.
# =============================================================================
set -uo pipefail
urls="$1"; bed="$2"; out="$3"; tmp="$4"; np="${5:-8}"
mkdir -p "$out" "$tmp"
export bed out tmp

one() {
    local acc="$1" url="$2" bw downloaded=0
    [ -s "$out/$acc.tsv.gz" ] && return 0
    if [[ "$url" =~ ^https?:// ]]; then
        bw="$tmp/$acc.bigWig"; downloaded=1
        if ! curl -sSf --retry 5 --retry-delay 20 -o "$bw.part" "$url"; then
            echo "$(date +%T) $acc FAILED download" >&2; rm -f "$bw.part"; return 1
        fi
        mv "$bw.part" "$bw"
    else
        bw="$url"                                   # local / mounted file: read in place, never removed
    fi
    if ! bigWigAverageOverBed "$bw" "$bed" "$tmp/$acc.tab"; then
        echo "$(date +%T) $acc FAILED bigWigAverageOverBed" >&2
        [ $downloaded -eq 1 ] && rm -f "$bw"
        rm -f "$tmp/$acc.tab"; return 1
    fi
    cut -f1,6 "$tmp/$acc.tab" | gzip -c > "$out/$acc.tsv.gz.part" && mv "$out/$acc.tsv.gz.part" "$out/$acc.tsv.gz"
    [ $downloaded -eq 1 ] && rm -f "$bw"
    rm -f "$tmp/$acc.tab"
    echo "$(date +%T) $acc done ($(ls "$out" | grep -c '\.tsv\.gz$') files)"
}
export -f one

cut -f1,2 "$urls" | xargs -P "$np" -n 2 bash -c 'one "$0" "$1"'
n_done=$(ls "$out" | grep -c '\.tsv\.gz$'); n_all=$(wc -l < "$urls")
echo "finished: $n_done of $n_all accessions"
[ "$n_done" -eq "$n_all" ]
