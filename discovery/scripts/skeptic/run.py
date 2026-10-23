"""Skeptic cycle 2: relax NaBiBr4 model vs Alexandria polymorphs / prototypes / decomposition products with
MACE-MPA-0 (+/- D3). Resumable: one json+cif per job in discovery/results/skeptic/relax/. Run from repo root."""
import json, sys, time, numpy as np
from pathlib import Path
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor as AA
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE
from mace.calculators import mace_mp

OUT = Path("discovery/results/skeptic/relax"); OUT.mkdir(parents=True, exist_ok=True)
ALEX = json.load(open("discovery/results/skeptic/alex.json"))
MODEL = "discovery/results/mace/REF-NaBiBr4-on-NaBiCl4proto/relaxed.cif"
_calc = {}
def calc(d3):
    if d3 not in _calc: _calc[d3] = mace_mp(model="medium-mpa-0", device="cpu", default_dtype="float64", dispersion=d3)
    return _calc[d3]

def alex(e):
    return Structure(Lattice(e["lattice_vectors"]), e["species_at_sites"], e["cartesian_site_positions"], coords_are_cartesian=True)

def job(name, s0, d3=False, const_vol=False, relax=True, steps=600):
    f = OUT/f"{name}.json"
    if f.exists(): return json.load(open(f))
    a = AA.get_atoms(s0); a.calc = calc(d3); n = len(a); t = time.time(); conv = None; ns = 0
    if relax:
        opt = FIRE(FrechetCellFilter(a, constant_volume=const_vol), logfile=None)
        conv = bool(opt.run(fmax=0.01, steps=steps)); ns = int(opt.nsteps)
        AA.get_structure(a).to(filename=str(OUT/f"{name}.cif"))
    r = dict(name=name, formula=s0.composition.reduced_formula, n=n, d3=d3, relaxed=relax, const_vol=const_vol,
             E=float(a.get_potential_energy()/n), V=float(a.get_volume()/n), V0=s0.volume/n, conv=conv, steps=ns, wall=time.time()-t)
    json.dump(r, open(f, "w")); print(r, flush=True); return r

SUB = {"Li": "Na", "K": "Na", "Rb": "Na", "Ag": "Na", "Sb": "Bi", "Cl": "Br", "I": "Br"}
def to_bromide(s):
    s = s.copy(); s.replace_species({k: v for k, v in SUB.items() if k in s.symbol_set})
    s.scale_lattice(33.0*len(s)); return s   # ~Alexandria NaBiBr4 volume/atom as a neutral start

if __name__ == "__main__":
    model = Structure.from_file(MODEL)
    for d3 in (False, True):
        t = "d3" if d3 else "nod"
        job(f"model_{t}", model, d3)
        for e in ALEX["BiBr4Na"] + ALEX["BiCl4Na"]:
            job(f"{e['id']}_{t}", alex(e), d3)
            if not d3: job(f"{e['id']}_sp", alex(e), relax=False)          # MACE at the PBE geometry
        for fam, ids in [("BrNa", ["agm003157587"]), ("BiBr3", ["agm002138960"]), ("BiBr6Na3", ["agm053446239"])]:
            for e in ALEX[fam]:
                if e["id"] in ids:
                    job(f"{e['id']}_{t}", alex(e), d3)
                    if not d3: job(f"{e['id']}_sp", alex(e), relax=False)
        # prototype enumeration: every neighbouring ABX4 entry within 30 meV of its hull, decorated as NaBiBr4
        for fam in ["BiCl4Na", "Br4NaSb", "Cl4NaSb", "BiI4Na", "I4NaSb", "BiBr4Li", "BiBr4K", "AgBiBr4", "BiBr4Rb"]:
            for e in ALEX[fam]:
                if e["_alexandria_hull_distance"] < 0.05 and e["nsites"] <= 24:
                    job(f"proto_{fam}_{e['id']}_{t}", to_bromide(alex(e)), d3)
    # fixed equal volumes (no D3): model vs the three Alexandria entries, shape+positions relaxed at constant V
    names = ["model"] + [e["id"] for e in ALEX["BiBr4Na"]]
    for v in (30.0, 33.0, 36.0, 39.0, 42.0):
        for nm in names:
            s = Structure.from_file(OUT/f"{nm}_nod.cif"); s.scale_lattice(v*len(s))
            job(f"fixV{v:.0f}_{nm}", s, False, const_vol=True, steps=400)
    print("DONE")
