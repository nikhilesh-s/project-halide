"""Strategy 3: local perturbations (supercell + rattle + isotropic strain) of the winner and of the KGaI4-/NaGaF4-type rivals,
looking for symmetry-lowering distortions. Cheap MACE for all; MACE+D3 directly for three winner supercells (s3_d3 jobs are
run by d3_final.py).  argv: minutes"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
SC = {"211": [2, 1, 1], "121": [1, 2, 1], "112": [1, 1, 2]}
def pert(parent, sc, sd, strain, seed):
    s = parent.copy()
    if sc: s.make_supercell(SC[sc])
    s.scale_lattice(s.volume*(1 + strain)**3); return rattle(s, sd, seed)
def plan():
    """(name, parent key, supercell, rattle A, linear strain, direct_d3)"""
    out = []
    for i, sc in enumerate(SC):
        out += [(f"pt_WIN_{sc}_r0.1_p3", "WIN_KBiF4type_DFTgeom", sc, 0.1, 0.03, False), (f"pt_WIN_{sc}_r0.3_m3", "WIN_KBiF4type_DFTgeom", sc, 0.3, -0.03, False),
                (f"pt_WIN_{sc}_r0.2_d3", "WIN_KBiF4type_DFTgeom", sc, 0.2, 0.03*(-1)**i, True)]
    for k in ("RIV_KGaI4type", "RIV_NaGaF4type"):
        t = k[4:-4]
        out += [(f"pt_{t}_111_r0.1_p3", k, "", 0.1, 0.03, False), (f"pt_{t}_111_r0.3_m3", k, "", 0.3, -0.03, False)]
        out += [(f"pt_{t}_{sc}_r0.2_{'pm'[i % 2]}3", k, sc, 0.2, 0.03*(-1)**i, False) for i, sc in enumerate(SC)]
    return out
def jobs(d3, dl):
    par = {k: S(json.load(open(C/"d3"/f"ref_{k}.json"))) for k in ("WIN_KBiF4type_DFTgeom", "RIV_KGaI4type", "RIV_NaGaF4type")}
    return [J(n, pert(par[k], sc, sd, st, 7 + i), d3, "s3_d3" if d3 else "s3_cheap", dl, fmax=0.01 if d3 else 0.03, steps=400)
            for i, (n, k, sc, sd, st, direct) in enumerate(plan()) if direct == d3]
if __name__ == "__main__":
    j = jobs(False, time.time() + float(sys.argv[1])*60)   # plan order: winner first, 48-atom rival supercells last
    print("S3 cheap finished", len(run(j, False, "s3_cheap", _relax2)), "of", len(j))
