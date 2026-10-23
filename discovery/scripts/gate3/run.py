"""Gate 3: saddle-point (rattle 0.03 A, seed 3) re-relaxation of cycle-3 hits and their Alexandria baselines/competitors.
MACE-MPA-0 float64 + D3(BJ), FIRE+FrechetCellFilter fmax 0.01. 3 single-thread workers, cells <= 16 atoms. Resumable.
Reuses discovery/results/skeptic3/relax/<name>.json when present (identical settings). Run from repo root."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import json, sys, time, itertools, urllib.request, urllib.parse
from pathlib import Path
from multiprocessing import Pool
OUT = Path("discovery/results/gate3/relax"); OUT.mkdir(parents=True, exist_ok=True)
S3 = Path("discovery/results/skeptic3/relax")
CHEMS = ["NaSbCl4", "NaBiCl4", "LiBiCl4", "LiSbBr4", "LiSbCl4", "LiBiBr4", "NaSbBr4", "LiAlBr4"]
MAXN = 16
U = "https://alexandria.icams.rub.de/pbe/v1/structures"
RF = "lattice_vectors,cartesian_site_positions,species_at_sites,nsites,chemical_formula_reduced,_alexandria_hull_distance,_alexandria_space_group"
def q(flt):
    url = U + "?" + urllib.parse.urlencode({"filter": flt, "page_limit": 100, "response_fields": RF}); out = []
    while url:
        for k in range(5):
            try: r = json.load(urllib.request.urlopen(url, timeout=90)); break
            except Exception as e: err = e; time.sleep(3*(k+1))
        else: raise err
        out += [{"id": d["id"], **d["attributes"]} for d in r["data"]]
        nx = (r.get("links") or {}).get("next"); url = nx.get("href") if isinstance(nx, dict) else nx; time.sleep(0.4)
    return out
def parse(c):
    import re
    A, M, X = re.match(r"(Li|Na)(Al|Sb|Bi)(Cl|Br)4", c).groups(); return A, M, X
def fetch():
    f = Path("discovery/results/gate3/alex.json")
    if f.exists(): return json.load(open(f))
    old = json.load(open("discovery/results/skeptic3/alex.json")); out = {}
    for c in CHEMS:
        if c in old: out[c] = old[c]; continue
        els = sorted(parse(c)); form = "".join(e + ("4" if e in ("Cl", "Br") else "") for e in els)
        abx4 = [e for e in q(f'chemical_formula_reduced="{form}"') if e["nsites"] <= 48]
        hull = {}
        for n in (1, 2, 3):
            for sub in itertools.combinations(els, n):
                has = ",".join(f'"{x}"' for x in sub)
                for e in q(f'elements HAS ALL {has} AND nelements={n} AND _alexandria_hull_distance<=0.005'): hull[e["id"]] = e
        out[c] = {"abx4": abx4, "hull": list(hull.values())}; print(c, len(abx4), len(hull), flush=True)
    json.dump(out, open(f, "w")); return out
_C = None
def job(spec):
    global _C
    name, sd, steps = spec
    for d in (OUT, S3):
        if (d/f"{name}.json").exists(): return json.load(open(d/f"{name}.json"))
    if _C is None:
        import torch; torch.set_num_threads(1)
        from mace.calculators import mace_mp
        _C = mace_mp(model="medium-mpa-0", device="cpu", default_dtype="float64", dispersion=True)
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor as AA
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    s0 = Structure.from_dict(sd); a = AA.get_atoms(s0); n = len(a); t = time.time(); a.calc = _C
    a.rattle(stdev=0.03, seed=3)
    opt = FIRE(FrechetCellFilter(a), logfile=None); conv = bool(opt.run(fmax=0.01, steps=steps))
    s1 = AA.get_structure(a); s1.to(filename=str(OUT/f"{name}.cif"))
    r = dict(name=name, formula=s0.composition.reduced_formula, n=n, E=float(a.get_potential_energy()/n), V=float(a.get_volume()/n),
             conv=conv, steps=int(opt.nsteps), wall=time.time()-t, structure=s1.as_dict())
    json.dump(r, open(OUT/f"{name}.json", "w")); print({k: v for k, v in r.items() if k != "structure"}, flush=True); return r
def alex(e):
    from pymatgen.core import Structure, Lattice
    return Structure(Lattice(e["lattice_vectors"]), e["species_at_sites"], e["cartesian_site_positions"], coords_are_cartesian=True)
def pick_hit(df, c):
    from pymatgen.core import Structure
    if c in ("NaSbCl4", "NaBiCl4", "LiBiCl4"):   # already DFT-tested / running: the cycle-3 TP-01 (rattled in skeptic3, reused by name)
        r = df[df.candidate == f"TP-{c}-01"].iloc[0]
        return r, Structure.from_dict(json.load(open(S3/f"TP-{c}-01_d3.json"))["structure"])
    h = df[(df.chem == c) & (df.status == "HIT_proposed_model") & (df.dE_D3_meV <= -5)].sort_values("E_decomp_meV")
    for _, r in h.iterrows():
        s = Structure.from_file(f"discovery/candidates/{r.candidate}/relaxed_mace_d3.cif"); p = s.get_primitive_structure()
        if len(p) <= 12: return r, p
    return h.iloc[0], None
def plan():
    import pandas as pd
    df = pd.read_csv("discovery/results/transplant_cycle3.csv"); A = fetch(); jobs = {}; sel = {}
    for c in CHEMS:
        r, p = pick_hit(df, c); sel[c] = dict(hit=r.candidate, best_existing_csv=r.best_existing_D3, natoms_prim=None if p is None else len(p)); 
        if p is None: continue
        jobs[f"{r.candidate}_d3"] = (f"{r.candidate}_d3", p.as_dict(), 600)
        ab = sorted(A[c]["abx4"], key=lambda e: e["_alexandria_hull_distance"]); keep = [e for e in ab if e["nsites"] <= MAXN][:3]
        keep += [e for e in ab if e["id"] == r.best_existing_D3 and e not in keep]
        for e in keep: jobs[f"{e['id']}_d3"] = (f"{e['id']}_d3", alex(e).as_dict(), 600)
        a_, m_, x_ = parse(c)
        for e in sorted(A[c]["hull"], key=lambda e: e["_alexandria_hull_distance"]):
            els = set(e["species_at_sites"]); nm = e["chemical_formula_reduced"]
            if e["nsites"] <= MAXN and (els in ({a_, x_}, {m_, x_}) or len(els) == 3):
                jobs.setdefault(f"{e['id']}_d3", (f"{e['id']}_d3", alex(e).as_dict(), 300))
        sel[c]["abx4"] = [(e["id"], e["nsites"], round(1000*e["_alexandria_hull_distance"], 1)) for e in ab]
    json.dump(sel, open("discovery/results/gate3/selection.json", "w"), indent=1); return jobs
if __name__ == "__main__":
    jobs = plan(); print("jobs", len(jobs), flush=True); t = time.time()
    with Pool(3) as p:
        for _ in p.imap_unordered(job, list(jobs.values())): pass
    print("DONE wall min", (time.time()-t)/60)
