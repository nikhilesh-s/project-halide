"""Skeptic cycle 3 relaxations. MACE-MPA-0 (medium-mpa-0, float64, CPU) +/- D3(BJ), FIRE + FrechetCellFilter, fmax 0.01,
no symmetry constraints -- same settings as discovery/scripts/skeptic/run.py. Every relaxation starts from a 0.03 A rattle
(seed 3) so baselines and hits get the same chance to leave a saddle. 4 single-thread workers. Resumable: one json+cif per
job in discovery/results/skeptic3/relax/. Run from repo root."""
import json, sys, time
from pathlib import Path
from multiprocessing import Pool
OUT = Path("discovery/results/skeptic3/relax"); OUT.mkdir(parents=True, exist_ok=True)
HITS = {"NaBiCl4": ["TP-NaBiCl4-01"], "NaSbCl4": ["TP-NaSbCl4-01", "TP-NaSbCl4-02"], "LiBiCl4": ["TP-LiBiCl4-01"]}
MAXN = 24
_C = {}
def calc(d3):
    if d3 not in _C:
        import torch; torch.set_num_threads(1)
        from mace.calculators import mace_mp
        _C[d3] = mace_mp(model="medium-mpa-0", device="cpu", default_dtype="float64", dispersion=d3)
    return _C[d3]

def job(spec):
    name, sd, d3, relax, steps = spec
    f = OUT/f"{name}.json"
    if f.exists(): return json.load(open(f))
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor as AA
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    s0 = Structure.from_dict(sd); a = AA.get_atoms(s0); n = len(a); t = time.time(); conv = None; ns = 0
    a.calc = calc(d3)
    if relax:
        a.rattle(stdev=0.03, seed=3)
        opt = FIRE(FrechetCellFilter(a), logfile=None)
        conv = bool(opt.run(fmax=0.01, steps=steps)); ns = int(opt.nsteps)
    s1 = AA.get_structure(a); s1.to(filename=str(OUT/f"{name}.cif"))
    r = dict(name=name, formula=s0.composition.reduced_formula, n=n, d3=d3, relaxed=relax, E=float(a.get_potential_energy()/n),
             V=float(a.get_volume()/n), V0=s0.volume/n, conv=conv, steps=ns, wall=time.time()-t, structure=s1.as_dict())
    json.dump(r, open(f, "w")); print({k: v for k, v in r.items() if k != "structure"}, flush=True); return r

def alex(e):
    from pymatgen.core import Structure, Lattice
    return Structure(Lattice(e["lattice_vectors"]), e["species_at_sites"], e["cartesian_site_positions"], coords_are_cartesian=True)

if __name__ == "__main__":
    from pymatgen.core import Structure
    A = json.load(open("discovery/results/skeptic3/alex.json")); jobs = {}; skipped = set()
    def add(name, s, d3, relax, steps=600): jobs.setdefault(name, (name, s.as_dict(), d3, relax, steps))
    for comp, hits in HITS.items():                       # hits + all existing ABX4 first (attacks 1, 2), then hull phases (3, 5)
        for h in hits:
            s = Structure.from_file(f"discovery/candidates/{h}/relaxed_mace_d3.cif")
            add(f"{h}_d3", s, True, True); add(f"{h}_nod", s, False, True)
        for e in A[comp]["abx4"]:
            if e["nsites"] > MAXN: skipped.add((e["id"], e["chemical_formula_reduced"], e["nsites"])); continue
            add(f"{e['id']}_d3", alex(e), True, True); add(f"{e['id']}_nod", alex(e), False, True); add(f"{e['id']}_sp", alex(e), False, False)
    for comp in HITS:
        for e in A[comp]["hull"] + A[comp].get("near", []):
            if e["nsites"] > MAXN: skipped.add((e["id"], e["chemical_formula_reduced"], e["nsites"])); continue
            add(f"{e['id']}_d3", alex(e), True, True, 300); add(f"{e['id']}_sp", alex(e), False, False)
    print("jobs", len(jobs), "skipped >24 atoms:", sorted(skipped), flush=True); t = time.time()
    with Pool(4) as p:
        for _ in p.imap_unordered(job, list(jobs.values())): pass
    print("DONE wall min", (time.time()-t)/60)
