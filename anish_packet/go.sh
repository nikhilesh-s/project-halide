#!/bin/bash
# Anish's one command. Safe to re-run any time (after the project owner says "updated", or once a morning):
#   bash anish_packet/go.sh          pull the latest PLAN, submit what is new in priority order, cancel what was dropped, print status
#   bash anish_packet/go.sh pack     end of window: sacct dump + one tarball of all results to send back
# Submits only the sbatch lines in PLAN (same commands as README). Never edits your code. Records each submission in
# $RUN_DIR/submitted/<tag> (Slurm job id), so nothing is submitted twice.
set -uo pipefail
cd "$(dirname "$0")"
git -C .. pull --no-rebase --no-edit -q || { echo "GIT PULL FAILED: send the project owner the error above and stop"; exit 1; }
source ./env.sh
sha256sum --quiet -c SHA256SUMS || { echo "PACKET CHECK FAILED: send the project owner the lines above and stop"; exit 1; }
cd "$RUN_DIR" && mkdir -p submitted logs

if [ "${1:-}" = pack ]; then
  sacct -u "$USER" -S 2026-10-16 --format=JobID,JobName%14,Partition,NodeList,AllocTRES%50,Elapsed,MaxRSS,State > sacct.txt
  tar -czf "$SCRATCH/anish_packet_results.tar.gz" -C "$(dirname "$RUN_DIR")" "$(basename "$RUN_DIR")"
  echo "PACKED: $SCRATCH/anish_packet_results.tar.gz ($(du -h "$SCRATCH/anish_packet_results.tar.gz" | cut -f1)). Send the project owner that file."; exit 0
fi

status() {
  echo "== queue"; squeue -h -u "$USER" -o "%.12i %.10j %.9T %.10M %R" | grep " pk_" || echo "   (no pk_ jobs queued or running)"
  echo "== DFT outputs finished (JOB DONE / inputs)"
  for d in "$PACKET"/S[0-9]/*/ "$PACKET"/S0/; do
    c=${d#$PACKET/}; c=${c%/}; n=$(ls "$d"*.in 2>/dev/null | wc -l); [ "$n" -gt 0 ] || continue
    k=$(grep -l "JOB DONE" "$RUN_DIR/$c"/*.out 2>/dev/null | wc -l); echo "   $c: $k/$n"
    if [ -f "$d/phonopy_disp.yaml" ] && [ "$k" -eq "$n" ] && [ ! -f "$RUN_DIR/$c/phonon_summary.json" ]; then
      "$ENV_BIN/python" "$PACKET/S1/post_phonons.py" "$PACKET/$c" "$RUN_DIR/$c" | tail -2
    fi
  done
}

# S0 control gate: nothing else is submitted until the cluster reproduces the laptop energy
S0OUT=$RUN_DIR/S0/scf_TP-NaBiCl4-01.out
if ! grep -q "JOB DONE" "$S0OUT" 2>/dev/null; then
  [ -s submitted/pk_s0 ] || sbatch --parsable --job-name=pk_s0 --ntasks=16 --mem=32G --time=00:30:00 --export=ALL,NK=4 \
    "$PACKET/run_pw.sh" S0 scf_TP-NaBiCl4-01 > submitted/pk_s0
  echo "S0 control submitted (job $(cat submitted/pk_s0)). Run go.sh again when it finishes (about 10 minutes)."; status; exit 0
fi
"$ENV_BIN/python" "$PACKET/S0/check_s0.py" "$S0OUT" || { echo "S0 FAILED: stop here and send the project owner the line above"; exit 1; }

while read -r kind tag after rest; do
  case "$kind" in ''|'#'*) continue ;; esac
  if [ "$kind" = cancel ]; then
    if [ -s "submitted/$tag" ] && [ ! -f "submitted/$tag.cancelled" ]; then
      scancel "$(cut -d';' -f1 "submitted/$tag")" && touch "submitted/$tag.cancelled" && echo "cancelled $tag (dropped from PLAN)"
    fi; continue
  fi
  if [ "$kind" = nice ]; then   # re-prioritise a set already submitted: "nice <tag> <value>"
    if [ -s "submitted/$tag" ] && [ "$(cat "submitted/$tag")" != manual ] && [ ! -f "submitted/$tag.nice$after" ]; then
      scontrol update JobId="$(cut -d';' -f1 "submitted/$tag")" Nice="$after" && touch "submitted/$tag.nice$after" && echo "re-prioritised $tag (nice $after)"
    fi; continue
  fi
  [ -s "submitted/$tag" ] && continue
  dep=""
  if [ "$after" != - ]; then
    [ -s "submitted/$after" ] || { echo "hold $tag: waits on $after"; continue; }
    dep="--dependency=afterok:$(cut -d';' -f1 "submitted/$after")"
  fi
  if [ "$kind" = ph ]; then
    read -r case time nice <<< "$rest"
    [ -n "$(squeue -h -u "$USER" -n "$tag")" ] && { echo "skip $tag: already in queue (submitted by hand)"; echo manual > "submitted/$tag"; continue; }
    n=$(ls "$PACKET/$case"/disp-*.in | wc -l)
    id=$(sbatch --parsable $dep --job-name="$tag" --array=0-$((n - 1))%6 --time="$time" --nice="$nice" "$PACKET/run_pw.sh" "$case")
  elif [ "$kind" = sb ]; then
    envs=${rest%%--*}; args=${rest#*-- }
    if { [[ $args =~ md_array\.sh\ ([^ ]+) ]] && [ -d "$RUN_DIR/${BASH_REMATCH[1]}" ]; } ||
       { [[ $args =~ run_pw\.sh\ ([^ ]+)\ ([^ ]+)$ ]] && [ -e "$RUN_DIR/${BASH_REMATCH[1]}/${BASH_REMATCH[2]}.out" ]; }; then
      echo "skip $tag: output already exists (submitted by hand)"; echo manual > "submitted/$tag"; continue
    fi
    id=$(eval "env $envs sbatch --parsable $dep $args")
  else
    echo "PLAN line not understood, skipped: $kind $tag"; continue
  fi
  [ -n "$id" ] && echo "$id" > "submitted/$tag" && echo "submitted $tag (job $id)"
done < "$PACKET/PLAN"
status
