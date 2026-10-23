#!/bin/bash
# Daytime runs on Arnav's M5, started from the dashboard: works through overnight/QUEUE in order with NP cores while the
# laptop stays usable. Finished jobs are skipped, relaxations resume from their last step, results are NOT pushed here
# (dashboard Push, Wi-Fi-off, or the overnight continue.sh do that).
#   NP=2 BG=1 bash overnight/m5_run.sh     light: 2 cores, background priority
#   NP=5 bash overnight/m5_run.sh          medium        NP=1 = offline-safe (no mpirun, no network needed)
cd "$(dirname "$0")/.."
if python3 -I -c "pass" >/dev/null 2>&1; then PY="python3 -I"; else PY="uv run --no-project python -I"; fi
[ "${BG:-0}" = 1 ] && taskpolicy -b -p $$   # background priority for this script and every job it starts
while true; do
  while read -r q; do
    [ -d "overnight/$q/jobs" ] && caffeinate -i $PY overnight/run_night.py "$q" --np "${NP:-2}" --stop "$(date -v-10M +%H:%M)" </dev/null >> overnight/m5_run.log 2>&1
  done < "${QUEUE:-overnight/QUEUE}"
  sleep 300   # retry failed launches (for example after the network changed)
done
