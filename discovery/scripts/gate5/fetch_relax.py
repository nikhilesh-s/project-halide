"""Gate 5 step 1: fetch the A3MX6 competing ternaries missing from the DFT queues from Alexandria (PBE) and give each
chosen cell the same MACE-MPA-0+D3 rattle/re-relax start as gates 3-4 (reuses gate4/run.py: q, alex, job; 1 thread).
Choice rule: lowest hull distance among cells <= 24 atoms; else the smallest cell <= 40 atoms (flagged); else none
("cannot be judged on laptop DFT"). ML relaxation only for cells <= 20 atoms; larger cells start from the Alexandria
PBE geometry. Run from repo root: taskpolicy -b discovery/.venv/bin/python -I discovery/scripts/gate5/fetch_relax.py"""
import importlib.util, json
from pathlib import Path
spec = importlib.util.spec_from_file_location("g4run", "discovery/scripts/gate4/run.py"); g4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(g4)
D = Path("discovery/results/gate5"); g4.OUT = D / "relax"; g4.OUT.mkdir(parents=True, exist_ok=True)
# competitor formula (Alexandria alphabetical reduced formula) -> candidate it judges
WANT = {"Cl6InLi3": "TP-LiInCl4-01", "Br6Li3Y": "TP-LiYBr4-02", "Br6Li3Sc": "TP-LiScBr4-01", "Br6Li3Sb": "TP-LiSbBr4-01",
        "BiBr6Li3": "TP-LiBiBr4-01", "BiCl6Li3": "TP-LiBiCl4-01", "Cl6Li3Sb": "TP-LiSbCl4-03", "Br6Na3Sb": "TP-NaSbBr4-01"}
f = D / "alex_competitors.json"
A = json.load(open(f)) if f.exists() else {}
for form in WANT:
    if form not in A: A[form] = g4.q(f'chemical_formula_reduced="{form}"'); json.dump(A, open(f, "w"))
sel = {}
for form, cand in WANT.items():
    ents = sorted(A[form], key=lambda e: e["_alexandria_hull_distance"])
    small = [e for e in ents if e["nsites"] <= 24]; mid = sorted((e for e in ents if e["nsites"] <= 40), key=lambda e: (e["nsites"], e["_alexandria_hull_distance"]))
    e = small[0] if small else (mid[0] if mid else None)
    sel[form] = dict(candidate=cand, listed=[(x["id"], x["nsites"], round(1000 * x["_alexandria_hull_distance"], 1), x["_alexandria_space_group"]) for x in ents],
                     chosen=e and e["id"], nsites=e and e["nsites"], hull_meV=e and round(1000 * e["_alexandria_hull_distance"], 1),
                     rule=None if e is None else ("lowest <=24 atoms" if small else "smallest cell <=40 atoms (lowest entries are larger)"))
    print(form, sel[form], flush=True)
    if e is None: continue
    s = g4.alex(e)
    if len(s.get_primitive_structure()) <= 20: g4.job((f"{e['id']}_d3", s.get_primitive_structure().as_dict(), 300))
    else: s.to(filename=str(g4.OUT / f"{e['id']}_d3.cif")); sel[form]["start"] = "Alexandria PBE geometry (no ML relax: > 20 atoms)"
json.dump(sel, open(D / "selection.json", "w"), indent=1)
# ternaries already in the gate-4 hull lists (<= 5 meV), for the record
G = json.load(open("discovery/results/gate4/alex.json"))
for c in ("LiInCl4", "LiScBr4", "LiYBr4", "NaGaI4"):
    print(c, "ternaries within 5 meV of Alexandria hull:", [(e["id"], e["chemical_formula_reduced"], e["nsites"], round(1000 * e["_alexandria_hull_distance"], 1)) for e in G[c]["hull"] if len(set(e["species_at_sites"])) == 3])
print("DONE")
