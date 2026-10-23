"""Attack 1: is the relaxed NaBiBr4 model the same structure as any Alexandria ABX4 entry? Run from repo root."""
import json, sys, numpy as np
from pathlib import Path
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
ALEX = json.load(open("discovery/results/skeptic/alex.json")); R = Path("discovery/results/skeptic/relax")
def alex(e): return Structure(Lattice(e["lattice_vectors"]), e["species_at_sites"], e["cartesian_site_positions"], coords_are_cartesian=True)
def red(s): return s.get_primitive_structure().get_reduced_structure("niggli")
SM = {"default": StructureMatcher(), "loose": StructureMatcher(ltol=0.3, stol=0.5, angle_tol=10),
      "vloose": StructureMatcher(ltol=0.4, stol=0.8, angle_tol=15)}
def cmp(a, b):
    out = []
    for k, m in SM.items():
        same = a.composition.reduced_formula == b.composition.reduced_formula
        ok = m.fit(a, b) if same else m.fit_anonymous(a, b)
        out.append(f"{k}={'Y' if ok else 'n'}")
    if same:
        d = SM["vloose"].get_rms_dist(a, b); out.append("rms=" + (f"{d[0]:.3f}" if d else "-"))
    return " ".join(out)
def desc(s, cut=1.25):
    """SG, and per-cation anion distances: CN = anions within cut x shortest cation-anion distance."""
    sg = [SpacegroupAnalyzer(s, symprec=p).get_space_group_symbol() for p in (0.01, 0.1)]
    o = {"sg": "/".join(sg), "V": round(s.volume/len(s), 2)}
    for el in sorted({str(x.specie) for x in s if str(x.specie) not in ("Cl", "Br", "I")}):
        seen = []
        for i, x in enumerate(s):
            if str(x.specie) != el: continue
            d = sorted(round(n.nn_distance, 2) for n in s.get_neighbors(x, 4.3) if str(n.specie) in ("Cl", "Br", "I"))
            cn = sum(v <= cut*d[0] for v in d); key = (cn, tuple(d[:8]))
            if key not in seen: seen.append(key)
        o[el] = [f"CN{c} {list(d)}" for c, d in seen]
    return o
model = red(Structure.from_file("discovery/results/mace/REF-NaBiBr4-on-NaBiCl4proto/relaxed.cif"))
print("MODEL", desc(model), "abc", [round(x, 2) for x in model.lattice.abc], [round(x, 1) for x in model.lattice.angles])
print("\n== vs Alexandria raw (PBE geometry) entries, all neighbouring ABX4 families ==")
for fam in ["BiBr4Na", "BiCl4Na", "Br4NaSb", "Cl4NaSb", "BiI4Na", "I4NaSb", "BiBr4Li", "BiBr4K", "AgBiBr4", "BiBr4Rb"]:
    for e in ALEX[fam]:
        s = red(alex(e)); print(fam, e["id"], f"hull={1000*e['_alexandria_hull_distance']:.1f}", f"sg#{e['_alexandria_space_group']}", f"n={len(s)}", cmp(model, s))
        if fam in ("BiBr4Na", "BiCl4Na"): print("    ", desc(s), "abc", [round(x, 2) for x in s.lattice.abc], [round(x, 1) for x in s.lattice.angles])
print("\n== vs MACE-relaxed NaBiBr4 structures (same settings) ==")
for tag in ("nod", "d3"):
    files = sorted(R.glob(f"*_{tag}.cif")); S = {}
    for f in files:
        s = Structure.from_file(f)
        if s.composition.reduced_formula == "NaBiBr4": S[f.stem] = red(s)
    m = S.get(f"model_{tag}")
    if m is None: continue
    print(tag, "model", desc(m))
    for k, s in S.items():
        if k != f"model_{tag}": print("  ", k, cmp(m, s), desc(s) if k.startswith("agm") else "")
    # dedupe all NaBiBr4 candidates into structure groups (loose)
    groups = SM["loose"].group_structures(list(S.values())); inv = {id(s): k for k, s in S.items()}
    print("  groups(loose):", [[inv[id(s)] for s in g] for g in groups])

def topo(s, cut=3.6):
    """Polyhedral connectivity: for each cation, sorted list of (neighbour cation, #shared anions, distance)."""
    an = ("Cl", "Br", "I"); c = cut*(s.volume/len(s)/38.9)**(1/3) if False else cut
    cat = [i for i, x in enumerate(s) if str(x.specie) not in an]; out = {}
    nb = {i: {(n.index, tuple(int(v) for v in n.image)) for n in s.get_neighbors(s[i], c) if str(n.specie) in an} for i in cat}
    for i in cat:
        L = []
        for n in s.get_neighbors(s[i], 8.0):
            if n.index not in cat: continue
            im = np.array(n.image, int)
            sh = len(nb[i] & {(j, tuple(int(v) for v in np.array(t)+im)) for j, t in nb[n.index]})
            if sh: L.append((str(n.specie), sh, round(n.nn_distance, 2)))
        out.setdefault(str(s[i].specie), set()).add(tuple(sorted(L)))
    return {k: [list(v) for v in vs] for k, vs in out.items()}
print("\n== connectivity (cation: neighbour cation, shared Br, distance A); Br cutoff 3.6 A ==")
for k in ["model_nod", "agm003733888_nod", "agm005739654_nod", "agm069058740_nod", "model_d3", "agm003733888_d3"]:
    f = R/f"{k}.cif"
    if f.exists(): print(k, json.dumps(topo(Structure.from_file(f))))
