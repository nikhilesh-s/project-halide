#!/bin/bash
#SBATCH --job-name=pk_pw
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
# One pw.x run, same launch line as the working dft/*/run*.sh jobs. Submit from $RUN_DIR after `source env.sh`.
#   single:  sbatch [overrides] $PACKET/run_pw.sh S0 scf_TP-NaBiCl4-01
#   array:   sbatch --array=1-18%6 $PACKET/run_pw.sh S1/TP-NaBiCl4-01      (task N runs disp-00N.in)
# NK (k-point pools) defaults to 2; use NK=1 if the job is killed for memory.
set -u
CASE=$1
NAME=${2:-$(printf "disp-%03d" "${SLURM_ARRAY_TASK_ID:?give an input name or use --array}")}
export PATH=$ENV_BIN:$PATH        # absolute env path: `conda activate` does not work inside Ionic batch jobs
export OMP_NUM_THREADS=1
OUT=$RUN_DIR/$CASE
QE_TMP=/tmp/${USER}_qe/$SLURM_JOB_ID        # node-local scratch; home quota is 5 GB
mkdir -p "$OUT" "$QE_TMP"
trap 'rm -rf "$QE_TMP"' EXIT
if grep -q "JOB DONE" "$OUT/$NAME.out" 2>/dev/null; then echo "$CASE/$NAME already done, skipping"; exit 0; fi
sed -e "s|__PSEUDO__|$PSEUDO_DIR|" -e "s|__TMP__|$QE_TMP|" "$PACKET/$CASE/$NAME.in" > "$OUT/$NAME.in"
echo "== $CASE/$NAME node=$(hostname) ntasks=$SLURM_NTASKS nk=${NK:-2} pw=$(which pw.x) $(date)"
mpirun --bind-to none -np "$SLURM_NTASKS" pw.x -nk "${NK:-2}" -i "$OUT/$NAME.in" > "$OUT/$NAME.out" 2>&1
grep -E "^!|Total force|JOB DONE|convergence NOT|Error|Estimated total dynamical RAM|PWSCF .*WALL" "$OUT/$NAME.out" | tail -n 6
echo "== end $(date)"
