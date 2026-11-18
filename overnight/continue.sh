#!/bin/bash
# One command for the M5: pull, make sure setup is current, run every queue in overnight/QUEUE in order
# until the stop time (default 07:15, 24h clock), then push results. Safe to re-run: finished jobs are skipped,
# unfinished relaxations resume from their last geometry.
#   bash overnight/continue.sh          # overnight, stops 07:15
#   bash overnight/continue.sh 18:00    # afternoon session, stops 18:00
set -uo pipefail
cd "$(dirname "$0")/.."
STOP="${1:-07:15}"
# some Macs block bare python3 with a shim; fall back to uv's python
if python3 -I -c "pass" >/dev/null 2>&1; then PY="python3 -I"; else PY="uv run --no-project python -I"; fi
git pull --no-rebase --no-edit || { echo "GIT PULL FAILED: fix that first, then re-run"; exit 1; }
bash overnight/setup.sh || { echo "SETUP FAILED: send the project owner the last 30 lines of overnight/build/make.log"; exit 1; }
PY="$PY" STOP="$STOP" caffeinate -i sh -c '
  while read -r q; do
    [ -d "overnight/$q/jobs" ] && $PY overnight/run_night.py "$q" --stop "$STOP" </dev/null  # mpirun reads stdin; without this it eats the rest of QUEUE
  done < overnight/QUEUE'
bash overnight/push.sh
