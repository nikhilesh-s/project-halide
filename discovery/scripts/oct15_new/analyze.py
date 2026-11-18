"""New-candidate (2026-10-15) analysis, copy of gate4/analyze.py: before/after-rattle dE and decomposition per chemistry (criteria as gate 3: dE <= -5 and
decomposition vs AX + MX3 <= -3 meV/atom after rattle). Also reports the hit against the MACE+D3 hull of every
relaxed competing phase (binaries and listed ternaries, hull <= 5 meV/atom in Alexandria). Run from repo root."""
import json, re, pandas as pd
from pathlib import Path
from pymatgen.core import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
R = {}
for d in ("skeptic3", "gate3", "gate4", "oct15_new"):
    for f in Path(f"discovery/results/{d}/relax").glob("*_d3.json"): R[f.stem[:-3]] = json.load(open(f))
A = json.load(open("discovery/results/oct15_new/alex.json")); sel = json.load(open("discovery/results/oct15_new/selection.json"))
df = pd.read_csv("discovery/results/transplant_cycle3.csv").set_index("candidate")
rows = []
for c, v in sel.items():
    A_, M_, X_ = re.match(r"(Li|Na)(Al|Ga|In|Sc|Y|Sb|Bi)(Cl|Br|I)4", c).groups(); h = v["hit"]; r0 = df.loc[h]
    amx4 = Composition(c).reduced_formula
    base = sorted((e["id"] for e in A[c]["abx4"] if e["id"] in R), key=lambda i: R[i]["E"]); bb = base[0]
    AX, MX, comp = {}, {}, {}
    for e in A[c]["hull"]:
        i = e["id"]
        if i not in R: continue
        f = Composition(R[i]["formula"]).reduced_formula
        if f == Composition(A_ + X_).reduced_formula: AX[i] = R[i]["E"]
        elif f == Composition(M_ + X_ + "3").reduced_formula: MX[i] = R[i]["E"]
        if f != amx4: comp[i] = f
    ax, mx = min(AX, key=AX.get), min(MX, key=MX.get); eh = R[h]["E"]
    # hull of all relaxed competitors; elements as dummies at 0 eV (all compound energies are negative totals)
    ents = [PDEntry(Composition(el), 0.0) for el in (A_, M_, X_)]
    ents += [PDEntry(Composition(f), R[i]["E"] * Composition(f).num_atoms, name=i) for i, f in comp.items()]
    pd_ = PhaseDiagram(ents); cc = Composition(c)
    dec_full = 1000 * (eh - pd_.get_hull_energy_per_atom(cc))
    route = {f"{e.name}({e.composition.reduced_formula})": round(x, 3) for e, x in pd_.get_decomposition(cc).items()}
    dE = 1000 * (eh - R[bb]["E"]); dec = 1000 * (eh - (2 * AX[ax] + 4 * MX[mx]) / 6)
    unrelaxed = [(e["id"], e["nsites"], round(1000 * e["_alexandria_hull_distance"], 1)) for e in A[c]["abx4"] if e["id"] not in R and e["_alexandria_hull_distance"] < 0.1]
    rows.append(dict(chem=c, hit=h, n=R[h]["n"], conv=R[h]["conv"], dE_before=r0.dE_D3_meV, dE_after=round(dE, 1), dec_before=r0.E_decomp_meV,
                     dec_after=round(dec, 1), dec_full_hull=round(dec_full, 1), route=route, best_after=bb, best_before=r0.best_existing_D3, AX=ax, MX3=mx,
                     baselines={i: (R[i]["n"], round(1000 * (R[i]["E"] - R[bb]["E"]), 1)) for i in base}, amx4_not_relaxed=unrelaxed,
                     competitors={i: (f, R[i]["n"]) for i, f in comp.items()}))
t = pd.DataFrame(rows); pd.set_option("display.width", 250, "display.max_colwidth", 120)
t["keep"] = (t.dE_after <= -5) & (t.dec_after <= -3)
print(t[["chem", "hit", "n", "conv", "dE_before", "dE_after", "dec_before", "dec_after", "dec_full_hull", "best_before", "best_after", "AX", "MX3", "keep"]].to_string())
for r in rows: print(r["chem"], "baselines", r["baselines"], "| not relaxed", r["amx4_not_relaxed"], "| full-hull route", r["route"])
t.to_json("discovery/results/oct15_new/summary.json", orient="records", indent=1)
