# Source this once per shell before any sbatch:   source anish_packet/env.sh
# Every path is overridable; defaults follow the layout seen in the cluster logs (scratch under /path/to/scratch/<user>).
export PACKET=${PACKET:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
export SCRATCH=${SCRATCH:-/path/to/scratch/$USER}
export ENV_BIN=${ENV_BIN:-$SCRATCH/matdisc_env/bin}          # env with pw.x (QE 7.5), mpirun, python, mace-torch, phonopy
export SSE_DIR=${SSE_DIR:-$HOME/matdiscovery/sse}            # folder holding the unmodified sse_md.py / sse_common.py
export PSEUDO_DIR=${PSEUDO_DIR:-$PACKET/pseudo}              # the exact UPF files of the laptop runs (md5 in pseudo/MD5SUMS)
export RUN_DIR=${RUN_DIR:-$SCRATCH/anish_packet_run}         # all outputs go here, never to the 5 GB home
export XDG_CACHE_HOME=${XDG_CACHE_HOME:-$SCRATCH/cache}      # MACE model cache
mkdir -p "$RUN_DIR/logs"
for x in "$ENV_BIN/pw.x" "$ENV_BIN/python" "$SSE_DIR/sse_md.py" "$PSEUDO_DIR/Bi.pbe-dn-kjpaw_psl.1.0.0.UPF"; do
  [ -e "$x" ] || echo "env.sh WARNING: missing $x"
done
