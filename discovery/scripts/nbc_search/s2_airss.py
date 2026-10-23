"""Strategy 2: AIRSS-style random structure search for NaBiCl4, Z = 1-4 (6-24 atoms in the primitive cell).
Random cell (26-30 A^3/atom, aspect <= 3:1, angles 60-120), random positions under minimum-distance constraints
(Na-Cl 2.4, Bi-Cl 2.3, Cl-Cl 3.0, cation-cation 3.2 A); half the starts are P1, half get a random triclinic / monoclinic /
orthorhombic space group (general positions, or snapped onto a symmetry element). Cheap MACE relax (float32, no D3), then
MACE+D3 float64 for the distinct lowest minima.  argv: n_starts cheap_minutes d3_minutes.  Resumable (seed = index)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from pymatgen.core import Lattice
from pymatgen.symmetry.groups import SpaceGroup
SGS = [2, 3, 4, 6, 7, 10, 11, 13, 14, 5, 8, 9, 12, 15] + list(range(16, 75))
CEN = {"P": 1, "C": 2, "A": 2, "I": 2, "F": 4}
_ops = {}
def ops(n):
    if n not in _ops:
        g = SpaceGroup.from_int_number(n); _ops[n] = ([(o.rotation_matrix, o.translation_vector) for o in g.symmetry_ops], CEN[g.symbol[0]], g.crystal_system)
    return _ops[n]
def dmin(a, b): return 3.0 if a == b == "Cl" else (2.4 if "Na" in (a, b) else 2.3) if "Cl" in (a, b) else 3.2

def gen(i):
    """Deterministic random start number i -> (Structure primitive, info) ."""
    rng = np.random.default_rng(1000 + i); scale = 1.0
    for attempt in range(2000):
        if attempt and attempt % 60 == 0: scale *= 0.95     # note: constraints relaxed 5 % per 60 failed attempts (dense Cl packing); recorded per start
        Z = int(rng.choice([1, 2, 3, 4], p=[0.15, 0.55, 0.05, 0.25])); n = 1 if rng.random() < 0.5 else int(rng.choice(SGS))
        O, cen, sysm = ops(n)
        if Z*cen*6 > 48: continue
        l = np.exp(rng.uniform(0, np.log(3), 3)); ang = [90, 90, 90]
        if sysm == "triclinic": ang = rng.uniform(60, 120, 3)
        elif sysm == "monoclinic": ang[1] = rng.uniform(90, 120)
        try: L = Lattice.from_parameters(*l, *ang); L = L.scale(rng.uniform(26, 30)*6*Z*cen)
        except Exception: continue
        if not np.isfinite(L.volume) or min(L.get_niggli_reduced_lattice().abc) < 3.3: continue
        fc, sp, fail = [], [], False
        for el, need in (("Bi", Z*cen), ("Na", Z*cen), ("Cl", 4*Z*cen)):
            tries = 0
            while need and not fail:
                tries += 1; fail = tries > 150
                p = rng.random(3)
                if len(O) > 1 and rng.random() < 0.5: R, t = O[rng.integers(1, len(O))]; p = (p + R @ p + t)/2   # onto a symmetry element when the op has a fixed point
                orb = []
                for R, t in O:
                    q = (R @ p + t) % 1
                    if not orb or L.get_all_distances(q, orb).min() > 0.5: orb.append(q)
                if len(orb) > need: continue
                if len(orb) > 1:
                    dd = L.get_all_distances(orb, orb) + 9*np.eye(len(orb))
                    if dd.min() < dmin(el, el)*scale: continue
                if fc:
                    dd = L.get_all_distances(orb, fc); lim = np.array([dmin(el, x) for x in sp])*scale
                    if (dd < lim).any(): continue
                fc += orb; sp += [el]*len(orb); need -= len(orb)
            if fail: break
        if fail: continue
        s = Structure(L, sp, fc).get_primitive_structure()
        if len(s) > 24 or len(s) % 6: continue
        return s, dict(Z=len(s)//6, sg_requested=n, constraint_scale=round(scale, 3), sg_start=sg(s, 0.01))
    raise RuntimeError("no start")

def distinct(res, tol=0.004):
    """res: name -> result dict. Clusters (lowest first): same primitive size, |dE| < tol, StructureMatcher match."""
    cl = []
    for n in sorted((n for n in res if ok(res[n])), key=lambda n: res[n]["E"]):
        p = prim(S(res[n])); e = res[n]["E"]
        c = next((c for c in cl if abs(c["E"]-e) < tol and same(p, c["p"])), None)
        if c: c["members"].append(n)
        else: cl.append(dict(rep=n, E=e, p=p, members=[n]))
    return cl

if __name__ == "__main__":
    N = int(sys.argv[1]); cm, dm = float(sys.argv[2]), float(sys.argv[3]); t = time.time()
    gf = C/"airss_starts.json"; G = json.load(open(gf)) if gf.exists() else {}
    for i in range(N):
        if f"rs_{i:04d}" not in G: s, info = gen(i); G[f"rs_{i:04d}"] = dict(info, structure=s.as_dict())
    json.dump(G, open(gf, "w")); print("starts", len(G), f"gen t={time.time()-t:.0f}s", flush=True)
    dl = time.time() + cm*60
    CH = run([J(n, Structure.from_dict(G[n]["structure"]), False, "airss_cheap", dl, steps=250) for n in sorted(G)[:N]], False, "airss_cheap", _relax2)
    good = {n: r for n, r in CH.items() if ok(r)}; print("cheap done", len(CH), "ok", len(good), f"t={time.time()-t:.0f}s", flush=True)
    print("S2 cheap finished")
