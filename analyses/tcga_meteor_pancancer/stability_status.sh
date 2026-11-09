#!/usr/bin/env bash
# Progress of the stability selection per cohort: chunks done / total, jobs in SLURM, summary written. Light: ./stability_status.sh
HERE=$(cd "$(dirname "$0")" && pwd); export PIPE_DIR=$HERE; source "$HERE/config/pipeline.conf"
SC=${METH_SCALE:-beta}
Q=$(command -v squeue >/dev/null && squeue -h -u "$USER" -o "%j %T" 2>/dev/null)
printf "%-7s %13s %-28s %s\n" cohort "chunks done" "SLURM jobs" summary
for C in $(cat "$WORK_DIR/cohorts_selected.txt" 2>/dev/null || ls "$WORK_DIR"); do
  D=$WORK_DIR/$C; [ -d "$D/chunks" ] || continue
  tot=$(ls $D/chunks/chunk_*.txt 2>/dev/null | wc -l); dn=$(ls $D/stability_$SC/done_chunk_* 2>/dev/null | wc -l)
  r=$(echo "$Q" | awk -v j="stab_$C" '$1==j && $2=="RUNNING"' | wc -l); pd=$(echo "$Q" | awk -v j="stab_$C" '$1==j && $2=="PENDING"' | wc -l)
  s=$(echo "$Q" | awk -v j="stabsum_$C" '$1==j{print $2}' | head -1)
  job="-"; [ $((r + pd)) -gt 0 ] && job="array: $r running, $pd pending"; [ -n "$s" ] && job="$job; summary $s"
  sm="-"; [ -f "$D/stability_$SC/stability_report_genes.tsv" ] && sm="written ($(($(wc -l < $D/stability_$SC/stability_report_genes.tsv) - 1)) genes)"
  [ "$dn" = 0 ] && [ "$job" = "-" ] && [ "$sm" = "-" ] && continue
  printf "%-7s %6s / %-5s %-28s %s\n" "$C" "$dn" "$tot" "$job" "$sm"
done
echo "logs: \$WORK_DIR/<COHORT>/logs/stability_*.out, stability_chunk_*.log; results: \$WORK_DIR/<COHORT>/stability_$SC/"
