#!/usr/bin/env bash
# Progress of every cohort: step, job state in SLURM, genes done / total, failures. Light, run anywhere:  ./status.sh
HERE=$(cd "$(dirname "$0")" && pwd); export PIPE_DIR=$HERE; source "$HERE/config/pipeline.conf"
SC=${METH_SCALE:-beta}
Q=$(command -v squeue >/dev/null && squeue -h -u "$USER" -o "%j %T" 2>/dev/null)
printf "%-7s %-10s %-30s %15s %8s %7s\n" cohort step "SLURM job" "genes done" "no data" failed
for C in $(cat "$WORK_DIR/cohorts_selected.txt" 2>/dev/null || ls "$WORK_DIR"); do
  D=$WORK_DIR/$C; [ -d "$D" ] || continue
  if   [ -f "$D/gene_list_$C.txt" ]; then step="genes"; tot=$(wc -l < "$D/gene_list_$C.txt")
  elif [ -f "$D/TCGA.$C.sampleMap/HumanMethylation450_$C" ]; then step="preprocessing"; tot="?"
  else step="not started"; tot="?"; fi
  prep=$(echo "$Q" | awk -v j="prep_$C" '$1==j{print $2}' | head -1)
  arr_r=$(echo "$Q" | awk -v j="meteor_$C" '$1==j && $2=="RUNNING"' | wc -l); arr_p=$(echo "$Q" | awk -v j="meteor_$C" '$1==j && $2=="PENDING"' | wc -l)
  done=$(ls "$D/results_$SC/performance" 2>/dev/null | wc -l)
  fail=$(cat "$D"/datasets_$SC/*/status.tsv 2>/dev/null | awk '$2=="failed"{print $1}' | sort -u | while read -r g; do [ -f "$D/results_$SC/performance/${g}_SVR.tsv" ] || echo "$g"; done | wc -l)   # still without a result
  nod=$(cat "$D"/datasets_$SC/*/status.tsv 2>/dev/null | awk '$2=="no_dataset"' | sort -u | wc -l)
  if [ -n "$prep" ]; then job="prep $prep"
  elif [ $((arr_r + arr_p)) -gt 0 ]; then job="genes: $arr_r running, $arr_p pending"
  elif [ "$step" = genes ] && [ $((done + nod)) -ge "$tot" ]; then job="-"; step="finished"
  elif [ "$step" = genes ] && [ $((done + nod + fail)) -ge "$tot" ]; then job="- (rerun the failed)"; step="finished"
  else job="NONE - stopped, see logs"; fi
  printf "%-7s %-10s %-30s %7s / %-6s %8s %7s\n" "$C" "$step" "$job" "$done" "$tot" "$nod" "$fail"
done
[ -n "$Q" ] && echo && echo "your jobs in the queue: $(echo "$Q" | grep -c .) ($(echo "$Q" | grep -c RUNNING) running)"
echo "no data = genes without a TSS anchor / probes in the window (not modelled, expected); failed genes are rerun by ./submit_cohort.sh <COHORT>"
echo "logs: $WORK_DIR/<COHORT>/logs/prepare_*.out (preparation), genes_*_*.out (gene tasks)"
