"""List every Alexandria PBE entry for given reduced formulas in a chemical system, sorted by hull distance."""
import json, sys, urllib.parse, urllib.request
B = "https://alexandria.icams.rub.de/pbe/v1/structures"
for cs, formula in [a.split(":") for a in sys.argv[1:]]:
    els = cs.split("-")
    flt = "elements HAS ALL " + ",".join(f'"{e}"' for e in els) + f" AND nelements={len(els)}"
    url = B + "?" + urllib.parse.urlencode({"filter": flt, "page_limit": 100,
        "response_fields": "chemical_formula_reduced,nsites,_alexandria_hull_distance,_alexandria_space_group,_alexandria_formation_energy_per_atom"})
    d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "scout"}), timeout=90))
    rows = sorted((a["_alexandria_hull_distance"], a["chemical_formula_reduced"], a["nsites"], a["_alexandria_space_group"], x["id"])
                  for x in d["data"] for a in [x["attributes"]] if formula in ("*", a["chemical_formula_reduced"]))
    print(f"== {cs} {formula}: {len(rows)} of {d['meta'].get('data_returned')} [{flt}]")
    for r in rows[:12]: print("  ehull=%.4f %s nsites=%s sg=%s %s" % r)
