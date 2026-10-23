"""Skeptic cycle 3 tables from discovery/results/skeptic3/relax/*.json -> stdout (saved as tables.txt). Run from repo root."""
import json, glob, sys, numpy as np
sys.path.insert(0, "discovery/scripts/skeptic3")
from pathlib import Path
from pymatgen.core import Structure, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from run import HITS, alex
R = {Path(f).stem: json.load(open(f)) for f in glob.glob("discovery/results/skeptic3/relax/*.json")}
A = json.load(open("discovery/results/skeptic3/alex.json"))
st = lambda k: Structure.from_dict(R[k]["structure"])
red = lambda s: s.get_primitive_structure().get_reduced_structure("niggli")
SM = {"default": StructureMatcher(), "loose": StructureMatcher(ltol=0.3, stol=0.5, angle_tol=10)}
AN = ("Cl",)
print("NOT CONVERGED:", [(k, r["formula"], r["steps"]) for k, r in R.items() if r["relaxed"] and not r["conv"]])

def desc(s, cut=1.25):   # as skeptic/identity.py: CN = anions within cut x shortest cation-anion distance
    o = {"sg": "/".join(SpacegroupAnalyzer(s, symprec=p).get_space_group_symbol() for p in (0.01, 0.1)), "V": round(s.volume/len(s), 2)}
    for el in sorted({str(x.specie) for x in s if str(x.specie) not in AN}):
        seen = []
        for x in s:
            if str(x.specie) != el: continue
            d = sorted(round(n.nn_distance, 2) for n in s.get_neighbors(x, 4.3) if str(n.specie) in AN)
            d = [float(v) for v in d]; key = (sum(v <= cut*d[0] for v in d), sum(v <= 3.3 for v in d), tuple(d[:9]))
            if key not in seen: seen.append(key)
        o[el] = [f"CN{c} (CN{c2} within 3.3 A) {list(d)}" for c, c2, d in seen]
    dm = s.distance_matrix + 99*np.eye(len(s)); sp = [str(x.specie) for x in s]
    o["shortest"] = {f"{a}-{b}": round(float(min(dm[i][j] for i in range(len(s)) for j in range(len(s)) if {sp[i], sp[j]} == {a, b} and (a != b or i != j) or (a == b and sp[i] == a and sp[j] == a))), 2)
                     for a, b in [("Cl", "Cl")] + [(c, "Cl") for c in sorted(set(sp)-{"Cl"})] + [tuple(sorted(set(sp)-{"Cl"}))]}
    return o

print("\n===== ATTACK 1: identity after rattle (0.03 A) + tight relaxation, MACE+D3 =====")
for comp, hits in HITS.items():
    ex = A[comp]["abx4"]
    for h in hits:
        r = R[f"{h}_d3"]; s = red(st(f"{h}_d3")); s0 = red(Structure.from_file(f"discovery/candidates/{h}/relaxed_mace_d3.cif"))
        print(f"\n{h}: conv={r['conv']} steps={r['steps']} E={r['E']:.5f} eV/atom; same as candidate CIF: default={SM['default'].fit(s, s0)} rms={SM['loose'].get_rms_dist(s, s0)[0]:.4f}")
        print("   ", json.dumps(desc(s)), "abc", [round(x, 2) for x in s.lattice.abc], [round(x, 1) for x in s.lattice.angles])
        for e in ex:
            row = []
            for tag, g in [("PBEgeom", red(alex(e)))] + ([("D3relaxed", red(st(f"{e['id']}_d3")))] if f"{e['id']}_d3" in R else []):
                d = SM["loose"].get_rms_dist(s, g); row.append(f"{tag}: default={'MATCH' if SM['default'].fit(s, g) else 'no'} loose={'MATCH' if SM['loose'].fit(s, g) else 'no'} rms={'%.3f' % d[0] if d else '-'}")
            print(f"    vs {e['id']} (sg#{e['_alexandria_space_group']}, n={e['nsites']}, hull +{1000*e['_alexandria_hull_distance']:.1f}): " + "; ".join(row))
    for a_, b_ in [(x, y) for i, x in enumerate(hits) for y in hits[i+1:]]:
        print(f"  {a_} vs {b_}: default={SM['default'].fit(red(st(a_+'_d3')), red(st(b_+'_d3')))} loose={SM['loose'].fit(red(st(a_+'_d3')), red(st(b_+'_d3')))}")
    print(f"  baseline descriptions ({comp}, D3 relaxed):")
    for e in ex: print("    ", e["id"], json.dumps(desc(red(st(f"{e['id']}_d3")))))

print("\n===== ATTACK 2: fair baseline (every Alexandria entry, rattled start, same settings) =====")
OLD = {}
for f in glob.glob("discovery/results/skeptic/relax/agm*_*.json") + glob.glob("discovery/results/transplant_cycle3/d3/ex_*.json"):
    j = json.load(open(f)); n = Path(f).stem.replace("ex_", ""); OLD[n if "_" in n else n+"_d3"] = j["E"]
for comp, hits in HITS.items():
    ids = [e["id"] for e in A[comp]["abx4"]]; hd = {e["id"]: 1000*e["_alexandria_hull_distance"] for e in A[comp]["abx4"]}
    b = {t: min(ids, key=lambda i: R[f"{i}_{t}"]["E"]) for t in ("d3", "nod")}
    print(f"\n{comp}: best existing D3 = {b['d3']}, best existing no-D3 = {b['nod']}; PBE hull distances {hd}")
    print("  structure | PBE hull (Alex) | MACE@PBEgeom rel | MACE relaxed rel (no D3) | MACE+D3 rel | V noD3 | V D3 | conv | dE(rattled - earlier unrattled run) D3/noD3 meV")
    rows = sorted(ids + hits, key=lambda i: R[f"{i}_d3"]["E"])
    for i in rows:
        g = lambda t: 1000*(R[f"{i}_{t}"]["E"]-R[f"{b[t if t != 'sp' else 'nod']}_{t}"]["E"])
        old = "/".join(f"{1000*(R[f'{i}_{t}']['E']-OLD[f'{i}_{t}']):+.2f}" if f"{i}_{t}" in OLD else "-" for t in ("d3", "nod"))
        print(f"  {i} | {('+%.1f' % hd[i]) if i in hd else '-'} | {('%+.1f' % g('sp')) if i in hd else '-'} | {g('nod'):+.1f} | {g('d3'):+.1f} | {R[f'{i}_nod']['V']:.2f} | {R[f'{i}_d3']['V']:.2f} | {R[f'{i}_d3']['conv']}/{R[f'{i}_nod']['conv']} | {old}")

print("\n===== ATTACK 3: MACE+D3 convex hull from Alexandria phases with PBE hull distance <= 5 meV/atom (+ all ABX4 entries) =====")
for comp, hits in HITS.items():
    ph = {e["id"]: e for e in A[comp]["hull"] + A[comp]["abx4"] + A[comp].get("near", []) if f"{e['id']}_d3" in R}
    ent = [PDEntry(st(f"{i}_d3").composition, R[f"{i}_d3"]["E"]*R[f"{i}_d3"]["n"], name=i) for i in ph]
    pd = PhaseDiagram(ent); c = Composition(comp)
    print(f"\n{comp}: {len(ent)} phases. Alexandria phases on the MACE+D3 hull:", sorted(f"{e.composition.reduced_formula}:{e.name}" for e in pd.stable_entries))
    print("   Alexandria phases, MACE+D3 e_above_hull (meV/atom) vs Alexandria PBE hull distance:")
    for i, e in sorted(ph.items(), key=lambda kv: kv[1]["chemical_formula_reduced"]):
        pe = next(x for x in ent if x.name == i); print(f"     {pe.composition.reduced_formula:9s} {i} n={e['nsites']:2d} PBE +{1000*e['_alexandria_hull_distance']:.1f} | MACE+D3 +{1000*pd.get_e_above_hull(pe):.1f} conv={R[i+'_d3']['conv']}")
    lo = lambda f: min((R[f"{i}_d3"]["E"] for i in ph if Composition(R[f"{i}_d3"]["formula"]).reduced_formula == Composition(f).reduced_formula), default=None)
    A_, M_ = comp[:2], comp[2:4]; ax, mx3, a3 = lo(A_+"Cl"), lo(M_+"Cl3"), lo(f"{A_}3{M_}Cl6")
    for h in hits + [min((e["id"] for e in A[comp]["abx4"]), key=lambda i: R[f"{i}_d3"]["E"])]:
        E = R[f"{h}_d3"]["E"]; dec, eh = pd.get_decomp_and_e_above_hull(PDEntry(c, E*c.num_atoms), allow_negative=True)
        s = f"   {h}: E - E_hull(Alexandria phases, MACE+D3) = {1000*eh:+.1f} meV/atom; hull products: " + " + ".join(f"{v:.3f} {k.composition.reduced_formula}[{k.name}]" for k, v in dec.items())
        s += f" | vs AX+MX3 {1000*(E-(2*ax+4*mx3)/6):+.1f}"
        if a3: s += f" | vs 1/3 A3MX6 + 2/3 MX3 {1000*(E-(10*a3+8*mx3)/18):+.1f}; (A3MX6 vs 3AX+MX3: {1000*(a3-(6*ax+4*mx3)/10):+.1f})"
        print(s)
        if h in hits:   # hull restricted to phases that keep M(III): AX, MX3, A3MX6, ABX4 (no redox with Cl2-rich phases / elements)
            keep = [x for x in ent if x.composition.reduced_formula in {Composition(f).reduced_formula for f in (A_+"Cl", M_+"Cl3", f"{A_}3{M_}Cl6", comp)}]
            from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram
            cpd = CompoundPhaseDiagram(keep + [PDEntry(c, E*c.num_atoms, name=h)], [Composition(A_+"Cl"), Composition(M_+"Cl3")])
            print("      pseudo-binary AX-MX3 hull incl. hit, vertices:", sorted(e.name for e in cpd.stable_entries))

print("\n===== ATTACK 5: model trust, MACE (no D3) vs Alexandria PBE =====")
for comp in HITS:
    ph = {e["id"]: e for e in A[comp]["hull"] + A[comp]["abx4"] + A[comp].get("near", []) if f"{e['id']}_sp" in R}
    for tag, sel in [("single point at PBE geometry, all phases", lambda i: f"{i}_sp"), ]:
        ent = {i: PDEntry(st(sel(i)).composition, R[sel(i)]["E"]*R[sel(i)]["n"], name=i) for i in ph}
        pd0 = PhaseDiagram([ent[i] for i, e in ph.items() if e["_alexandria_hull_distance"] == 0.0])      # hull of Alexandria's own hull-0 phases, MACE energies
        rows = []
        for i, e in ph.items():
            m = 1000*pd0.get_decomp_and_e_above_hull(ent[i], allow_negative=True, check_stable=False)[1] if e["_alexandria_hull_distance"] > 0 else None
            if m is not None: rows.append((e["chemical_formula_reduced"], i, 1000*e["_alexandria_hull_distance"], m))
        err = np.array([m-p for _, _, p, m in rows]); t = [r for r in rows if len(Composition(r[0]).elements) == 3]; et = np.array([m-p for _, _, p, m in t])
        print(f"\n{comp} [{tag}]: {len(rows)} non-hull entries; RMS(MACE - PBE hull distance) = {np.sqrt((err**2).mean()):.1f} meV/atom, max |err| {np.abs(err).max():.1f}; ternaries only (n={len(t)}): RMS {np.sqrt((et**2).mean()):.1f}, max {np.abs(et).max():.1f}")
        for f, i, p, m in sorted(rows): print(f"     {f:9s} {i}  PBE +{p:.1f}  MACE {m:+.1f}  err {m-p:+.1f}")
        neg = [(ent[i].composition.reduced_formula, i, round(float(1000*PhaseDiagram([v for k, v in ent.items() if k != i]).get_decomp_and_e_above_hull(ent[i], allow_negative=True)[1]), 1)) for i, e in ph.items() if e["_alexandria_hull_distance"] == 0.0 and len(ent[i].composition.elements) > 1]
        print("     Alexandria hull-0 compounds, MACE energy vs MACE hull of the others (should be <= 0):", neg)
    ids = [e["id"] for e in A[comp]["abx4"]]; hd = {e["id"]: 1000*e["_alexandria_hull_distance"] for e in A[comp]["abx4"]}
    for t in ("sp", "nod"):
        d = [(1000*(R[f"{i}_{t}"]["E"]-R[f"{j}_{t}"]["E"]))-(hd[i]-hd[j]) for k, i in enumerate(ids) for j in ids[k+1:]]
        print(f"   ABX4 polymorph pair differences ({'PBE geometry' if t == 'sp' else 'MACE relaxed'}), n={len(d)} pairs: RMS err {np.sqrt(np.mean(np.square(d))):.1f}, max {np.abs(d).max():.1f} meV/atom")

print("\n===== ATTACK 4: phonons (2x2x2, commensurate q only, MACE+D3) =====")
for f in sorted(glob.glob("discovery/results/skeptic3/phonon/*.json")):
    p = json.load(open(f)); s = f"{p['id']}: n_disp={p['n_disp']} min={p['min_THz']:.2f} THz at q={p['q_min']}; modes < -0.3 THz: {p['n_below_m0p3']}; per-q min: " + " ".join(f"{v:.2f}" for v in p["per_q_min"].values())
    for k, v in p.get("follow", {}).items():
        g = st(f"{p['id']}_follow{k}_d3"); s += f"\n     follow {k}: n={v['n']} dE={v['dE_meV']:+.2f} meV/atom conv={v['conv']} SG {desc(red(g))['sg']} V={v['V']:.2f}; matches parent (default): {SM['default'].fit(red(g), red(st(p['id']+'_d3')))}"
    print(s)
E = lambda k: R[k+"_d3"]["E"]
print("\nRedox check, MACE+D3: 5 NaSbCl4 -> 3 NaSbCl6 + 2 Sb + 2 NaCl (30 atoms). Product energy minus reactant (meV/atom of NaSbCl4), positive = NaSbCl4 favoured:")
prod = (3*8*E("agm003219938")+2*E("agm003157733")+2*2*E("agm003157609"))/30
for h in ("TP-NaSbCl4-01", "TP-NaSbCl4-02", "agm005743441"): print(f"   {h}: {1000*(prod-E(h)):+.1f}")
