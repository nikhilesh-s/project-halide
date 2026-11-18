"""Gate 3 analysis: before/after-rattle dE and decomposition per chemistry. Run from repo root."""
import json, re, pandas as pd
from pathlib import Path
R = {}
for d in (Path("discovery/results/skeptic3/relax"), Path("discovery/results/gate3/relax")):
    for f in d.glob("*_d3.json"): R[f.stem[:-3]] = json.load(open(f))
A = json.load(open("discovery/results/gate3/alex.json")); sel = json.load(open("discovery/results/gate3/selection.json"))
df = pd.read_csv("discovery/results/transplant_cycle3.csv").set_index("candidate")
rows = []; detail = {}
for c, v in sel.items():
    A_, M_, X_ = re.match(r"(Li|Na)(Al|Sb|Bi)(Cl|Br)4", c).groups(); h = v["hit"]; r0 = df.loc[h]
    ab = [e["id"] for e in A[c]["abx4"] if e["id"] in R]
    base = sorted(ab, key=lambda i: R[i]["E"]); bb = base[0]
    # state which ABX4 baselines were relaxed
    ent = {e["id"]: e for e in A[c]["hull"]}; AX, MX, TER = {}, {}, {}
    for i, e in ent.items():
        if i not in R: continue
        els = set(e["species_at_sites"]); f = e["chemical_formula_reduced"]
        if els == {A_, X_}: AX[i] = R[i]["E"]
        elif els == {M_, X_}: MX[i] = R[i]["E"]
        elif len(els) == 3 and f != R[bb]["formula"]: TER[i] = (f, R[i]["E"], e["nsites"], round(1000*e["_alexandria_hull_distance"], 1))
    ax, mx = min(AX, key=AX.get), min(MX, key=MX.get); eh = R[h]["E"]
    dE = 1000*(eh - R[bb]["E"]); dec = 1000*(eh - (2*AX[ax] + 4*MX[mx])/6)
    t3 = {i: t for i, t in TER.items() if re.fullmatch(rf"{A_}3{M_}{X_}6|{X_}6{A_}3{M_}|{A_}3{M_}{X_}6", t[0]) or t[0].startswith(A_+"3")}
    dect = {i: 1000*(eh - (10/3*t[1] + 8/3*MX[mx])/6) for i, t in t3.items()} if t3 else {}
    rows.append(dict(chem=c, hit=h, n=R[h]["n"], dE_before=r0.dE_D3_meV, dE_after=round(dE, 1), dec_before=r0.E_decomp_meV, dec_after=round(dec, 1),
                     best_after=bb, best_before=r0.best_existing_D3, AX=ax, MX3=mx, dec_vs_A3MX6={k: round(x, 1) for k, x in dect.items()},
                     conv=R[h]["conv"], baselines={i: round(1000*(R[i]["E"] - R[bb]["E"]), 1) for i in base}))
    detail[c] = dict(ternaries=TER)
t = pd.DataFrame(rows); pd.set_option("display.width", 250, "display.max_colwidth", 80)
t["keep"] = (t.dE_after <= -5) & (t.dec_after <= -3)
print(t.drop(columns=["baselines"]).to_string()); print(t[["chem", "baselines"]].to_string())
for c, d in detail.items(): print(c, d["ternaries"])
t.to_json("discovery/results/gate3/summary.json", orient="records", indent=1)
