#!/bin/sh
# Runs the stages one after another (never more than 4 workers). Usage: sh chain.sh  (from repo root; everything is resumable)
PY="taskpolicy -b discovery/.venv/bin/python -I -W ignore discovery/scripts/nbc_search"
L=discovery/cache/nbc_search
while pgrep -f "nbc_search/s2_airss.py" >/dev/null || pgrep -f "s1_fetch.py" >/dev/null; do sleep 5; done
$PY/s1_run.py 13 > $L/s1.log 2>&1
$PY/s3_perturb.py 5 > $L/s3.log 2>&1
$PY/d3_final.py 15 > $L/d3.log 2>&1
$PY/phonon.py 6 > $L/ph.log 2>&1
$PY/report.py > $L/report.log 2>&1
