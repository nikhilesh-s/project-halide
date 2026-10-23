#!/bin/bash
#SBATCH --job-name=pk_md
#SBATCH --partition=all
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
# Same layout as code/md_array.sh (NPROC python workers share one GPU), but calls the trajectory-saving wrapper.
# Submit from $RUN_DIR after `source env.sh`:
#   NTASK=<n> NPROC=4 sbatch --array=0-<n-1> $PACKET/S2/md_array.sh <subdir> --items <items.json.gz> --tag <tag>
# n = ceil(number of items / NPROC). Results: $RUN_DIR/<subdir>/md_<tag>_mace_<task>_<proc>.json.gz and traj/*.npz
set -u
export PATH=$ENV_BIN:$PATH
export OMP_NUM_THREADS=2
SUB=$1; shift
mkdir -p "$RUN_DIR/$SUB/logs" && cd "$RUN_DIR/$SUB"
export TRAJ_DIR=$RUN_DIR/$SUB/traj
NPROC=${NPROC:-4}
echo "node=$(hostname) gpu=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader) task=$SLURM_ARRAY_TASK_ID python=$(which python)"
for p in $(seq 0 $((NPROC-1))); do
  python "$PACKET/S2/md_traj_wrapper.py" --proc $p --nproc $NPROC "$@" > logs/md_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_p${p}.log 2>&1 &
done
( while sleep 60; do nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv,noheader; done > logs/vram_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.csv ) &
VP=$!
wait $(jobs -p | grep -v "^$VP$")
kill $VP
echo "all procs done"
