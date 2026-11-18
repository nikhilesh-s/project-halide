"""Skeptic cycle 3: fetch Alexandria PBE entries. (a) every ABX4 entry of the three hit compositions (any hull, <=48 atoms);
(b) every phase with hull distance <= 5 meV/atom in each ternary chemical system (elements, binaries, ternaries).
Cache: discovery/results/skeptic3/alex.json. Run from repo root."""
import json, time, urllib.request, urllib.parse, itertools
U = "https://alexandria.icams.rub.de/pbe/v1/structures"
SYS = {"NaBiCl4": ("Na", "Bi", "Cl"), "NaSbCl4": ("Na", "Sb", "Cl"), "LiBiCl4": ("Li", "Bi", "Cl")}
RF = "lattice_vectors,cartesian_site_positions,species_at_sites,nsites,chemical_formula_reduced,_alexandria_hull_distance,_alexandria_space_group,_alexandria_energy,_alexandria_formation_energy_per_atom"
def q(flt):
    url = U + "?" + urllib.parse.urlencode({"filter": flt, "page_limit": 100, "response_fields": RF}); out = []
    while url:
        for k in range(5):
            try: r = json.load(urllib.request.urlopen(url, timeout=90)); break
            except Exception as e: err = e; time.sleep(3*(k+1))
        else: raise err
        out += [{"id": d["id"], **d["attributes"]} for d in r["data"]]
        nx = (r.get("links") or {}).get("next"); url = nx.get("href") if isinstance(nx, dict) else nx; time.sleep(0.5)
    return out
rform = lambda els: "".join(e + ("4" if e == "Cl" else "") for e in sorted(els))
out = {}
for name, els in SYS.items():
    abx4 = [e for e in q(f'chemical_formula_reduced="{rform(els)}"') if e["nsites"] <= 48]
    hull = {}
    for n in (1, 2, 3):
        for sub in itertools.combinations(els, n):
            has = ",".join(f'"{x}"' for x in sub)
            for e in q(f'elements HAS ALL {has} AND nelements={n} AND _alexandria_hull_distance<=0.005'): hull[e["id"]] = e
    has = ",".join(f'"{x}"' for x in els)   # (c) extended competitors: every ternary within 50 meV/atom of the PBE hull
    near = [e for e in q(f'elements HAS ALL {has} AND nelements=3 AND _alexandria_hull_distance<=0.05') if e["id"] not in hull and e["id"] not in {x["id"] for x in abx4}]
    out[name] = {"abx4": abx4, "hull": list(hull.values()), "near": near}
    print("   other ternaries <=50meV:", sorted((e["chemical_formula_reduced"], e["id"], round(1000*e["_alexandria_hull_distance"], 1), e["nsites"]) for e in near))
    print(name, "ABX4:", [(e["id"], e["_alexandria_space_group"], round(1000*e["_alexandria_hull_distance"], 1), e["nsites"]) for e in abx4])
    print("   hull<=5meV:", sorted((e["chemical_formula_reduced"], e["id"], e["_alexandria_space_group"], round(1000*e["_alexandria_hull_distance"], 1), e["nsites"]) for e in hull.values()))
json.dump(out, open("discovery/results/skeptic3/alex.json", "w"))
