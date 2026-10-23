"""Cycle-4 scout: Materials Cloud MC3D (structures sourced from ICSD/COD/MPDS) OPTIMADE lookup per ternary chemsys."""
import json, sys, urllib.request, urllib.parse

for cs in sys.argv[1:]:
    f = "elements HAS ALL " + ",".join('"%s"' % e for e in cs.split("-")) + " AND nelements=3"
    u = "https://optimade.materialscloud.org/main/mc3d-pbe-v1/v1/structures?" + urllib.parse.urlencode({"filter": f, "page_limit": 50})
    try:
        d = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "curl/8"}), timeout=40))["data"]
        print("MC3D", cs, len(d))
        for x in d:
            a = x["attributes"]
            extra = {k: v for k, v in a.items() if k.startswith("_") and v is not None and not isinstance(v, (list, dict))}
            print("  ", x["id"], a.get("chemical_formula_reduced"), "nsites=%s" % a.get("nsites"), extra)
    except Exception as e:
        print("MC3D", cs, "ERR", type(e).__name__, str(e)[:100])
