"""Shared bits for the NaBiCl4 adversarial search (cycle 5). Reuses the relaxation worker of transplant_scan.py
(MACE-MPA-0 medium, FIRE + FrechetCellFilter; cheap = float32 no D3, d3 = float64 + D3(BJ)). Everything is a PROPOSED MODEL.
Run every script from the repo root with: taskpolicy -b discovery/.venv/bin/python -I discovery/scripts/nbc_search/<script>.py"""
import json, sys, time, os
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, "discovery/scripts"); sys.path.insert(0, "overnight")
from transplant_scan import _init, _relax, struct, prim, same, ok, SM   # noqa: reuse, do not re-implement
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

C = Path("discovery/cache/nbc_search"); C.mkdir(parents=True, exist_ok=True)
NPROC = 4
DFT_OUT = "overnight/night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.out"
ALEX5 = ["agm005743614", "agm088474086", "agm084661913", "agm003736714", "agm075995748"]

def rattle(s, sd, seed):
    s = s.copy(); rng = np.random.default_rng(seed)
    for i, d in enumerate(rng.normal(0, sd, (len(s), 3))): s.translate_sites([i], d, frac_coords=False)
    return s

def J(name, s, d3, tag, deadline, fmax=None, steps=None):
    return (name, s.as_dict(), fmax or (0.01 if d3 else 0.05), steps or (600 if d3 else 300), str(C/tag), deadline)

def _relax2(job):
    """Same contract as transplant_scan._relax but LBFGS instead of FIRE (fewer steps; used for the cheap triage stages only,
    because `taskpolicy -b` makes every step ~10x slower than in the earlier scans)."""
    import transplant_scan as T
    name, sd, fmax, steps, out, deadline = job
    f = Path(out)/f"{name}.json"
    if f.exists(): return json.load(open(f))
    if time.time() > deadline: return None
    from pymatgen.io.ase import AseAtomsAdaptor as AA
    from ase.filters import FrechetCellFilter
    from ase.optimize import LBFGS
    s0 = Structure.from_dict(sd); a = AA.get_atoms(s0); a.calc = T._C; n = len(a); t = time.time()
    try:
        opt = LBFGS(FrechetCellFilter(a), logfile=None, maxstep=0.15)
        conv = bool(opt.run(fmax=fmax, steps=steps))
        s1 = AA.get_structure(a); dm = s1.distance_matrix + np.eye(n)*1e6
        r = dict(name=name, n=n, E=float(a.get_potential_energy()/n), V=s1.volume/n, V0=s0.volume/n, conv=conv,
                 steps=int(opt.nsteps), min_dist=float(dm.min()), wall=time.time()-t, structure=s1.as_dict())
    except Exception as e:
        r = dict(name=name, n=n, error=repr(e)[:200], wall=time.time()-t)
    json.dump(r, open(f, "w")); return r

def run(jobs, d3, tag, fn=_relax):
    """jobs from J(); resumable (one json per job); jobs not started before their deadline are skipped."""
    out = C/tag; out.mkdir(exist_ok=True); t = time.time()
    todo = [j for j in jobs if not (out/f"{j[0]}.json").exists()]
    print(f"[{tag}] {len(jobs)} jobs, {len(todo)} to run", flush=True)
    if todo:
        with Pool(NPROC, initializer=_init, initargs=(d3,)) as p:
            for i, r in enumerate(p.imap_unordered(fn, todo, chunksize=1)):
                if i % 25 == 0: print(f"  [{tag}] {i}/{len(todo)} t={time.time()-t:.0f}s", flush=True)
    return {j[0]: json.load(open(out/f"{j[0]}.json")) for j in jobs if (out/f"{j[0]}.json").exists()}

def S(r): return Structure.from_dict(r["structure"])
def sg(s, sp=0.1):
    try: return SpacegroupAnalyzer(s, symprec=sp).get_space_group_symbol()
    except Exception: return "?"

def refs():
    """name -> starting structure of the nine reference NaBiCl4 structures."""
    import make_inputs as mi
    A = {e["id"]: struct(e) for e in json.load(open("discovery/cache/alex_BiCl4Na.json"))}
    R = {"WIN_KBiF4type_DFTgeom": mi.final_structure(DFT_OUT),
         "WIN_KBiF4type_MLgeom": Structure.from_file("discovery/results/skeptic3/relax/TP-NaBiCl4-01_d3.cif"),
         "RIV_KGaI4type": Structure.from_file("discovery/candidates/TP-NaBiCl4-01/relaxed_mace_d3.cif"),
         "RIV_NaGaF4type": Structure.from_file("discovery/candidates/TP-NaBiCl4-03/relaxed_mace_d3.cif")}
    R.update({i: A[i] for i in ALEX5})
    R["agm084661913_sk3geom"] = Structure.from_file("discovery/results/skeptic3/relax/agm084661913_d3.cif")   # its rattle-found Fdd2 minimum
    return R

def ref_results():
    """Relaxed refs: (cheap dict, d3 dict). Rattled 0.03 A first (lesson of skeptic cycle 3)."""
    R = refs(); dl = time.time() + 1800
    d3 = run([J("ref_" + k, rattle(s, 0.03, 3), True, "d3", dl) for k, s in R.items()], True, "d3")
    ch = run([J("ref_" + k, rattle(s, 0.03, 3), False, "cheap", dl, fmax=0.02, steps=500) for k, s in R.items()], False, "cheap")
    return {k[4:]: v for k, v in ch.items()}, {k[4:]: v for k, v in d3.items()}

def identify(s, lib):
    """lib: name -> list of structures. Returns 'name(default|loose)' of the first match or ''."""
    p = prim(s)
    for n, ss in lib.items():
        for x in ss:
            if (t := same(p, x)): return f"{n}({t})"
    return ""

def ref_lib(ch, d3):
    R = refs(); return {k: [prim(R[k])] + [prim(S(d[k])) for d in (ch, d3) if k in d and ok(d[k])] for k in R}

__all__ = [k for k in dir() if not k.startswith("__")]   # scripts do `from common import *`, underscore names included

if __name__ == "__main__":   # self-check + calibration print
    t = time.time(); ch, d3 = ref_results(); w = d3["WIN_KBiF4type_DFTgeom"]["E"]
    for k, v in d3.items(): print(f"{k:26s} n={v['n']:2d} D3 {1000*(v['E']-w):+7.2f}  V {v['V']:.2f} conv {v['conv']} steps {v['steps']} wall {v['wall']:.0f}s | cheap {1000*(ch[k]['E']-ch['WIN_KBiF4type_DFTgeom']['E']):+7.2f} V {ch[k]['V']:.2f} steps {ch[k]['steps']} wall {ch[k]['wall']:.0f}s  sg {sg(S(v))}")
    L = ref_lib(ch, d3); assert identify(S(d3["WIN_KBiF4type_MLgeom"]), {"w": L["WIN_KBiF4type_DFTgeom"]}), "winner ML/DFT geometries must match"
    print("calibration: MACE+D3 winner - agm088474086 = %.2f meV/atom (DFT -5.7)" % (1000*(w - d3["agm088474086"]["E"])), "wall", time.time()-t)
