"""Strategy 1, step 2: cheap MACE pass (float32, no D3, fmax 0.05, <= 150 steps) on the new prototypes transplanted to NaBiCl4
at 28 A^3/atom, <= 24 atoms first, then 20 meV hull bins, then smallest cells, until the deadline.  argv: minutes"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
if __name__ == "__main__":
    P = sorted(json.load(open(C/"protos.json"))["protos"], key=lambda p: (p["n"] > 24, int(p["hull"]/0.02), p["n"], p["hull"]))
    dl = time.time() + float(sys.argv[1])*60
    R = run([J("tp_" + p["id"], Structure.from_dict(p["structure"]), False, "s1_cheap", dl, steps=150) for p in P], False, "s1_cheap", _relax2)
    print("S1 cheap finished", len(R), "of", len(P))
