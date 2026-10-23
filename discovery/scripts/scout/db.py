import json, os, urllib.request, urllib.parse

key = open(os.path.expanduser("~/.mp_api_key")).read().strip()
for cs in ["Na-Bi-Cl", "Na-Sb-Cl", "Li-Bi-Cl"]:
    u = "https://api.materialsproject.org/materials/summary/?" + urllib.parse.urlencode(
        {"chemsys": cs, "_fields": "material_id,formula_pretty,energy_above_hull,theoretical,symmetry", "_limit": 200})
    try:
        r = urllib.request.Request(u, headers={"X-API-KEY": key, "accept": "application/json", "User-Agent": "curl/8"})
        d = json.load(urllib.request.urlopen(r, timeout=60))["data"]
        print("MP", cs, len(d))
        for x in sorted(d, key=lambda x: x["formula_pretty"]):
            print("  ", x["material_id"], x["formula_pretty"], round(x["energy_above_hull"], 4), x["theoretical"], x["symmetry"]["symbol"])
    except Exception as e:
        print("MP", cs, "ERR", type(e).__name__, str(e)[:100])

for els, label in [(("Na", "Bi", "Cl"), "Na-Bi-Cl"), (("Na", "Sb", "Cl"), "Na-Sb-Cl"), (("Li", "Bi", "Cl"), "Li-Bi-Cl")]:
    f = 'elements HAS ALL ' + ",".join('"%s"' % e for e in els) + ' AND nelements=3'
    u = "https://oqmd.org/optimade/structures?" + urllib.parse.urlencode(
        {"filter": f, "page_limit": 50, "response_fields": "chemical_formula_reduced,_oqmd_entry_id,_oqmd_spacegroup,_oqmd_stability"})
    try:
        d = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "curl/8"}), timeout=180))
        print("OQMD", label, len(d.get("data", [])))
        for x in d.get("data", []):
            a = x["attributes"]
            print("  ", x["id"], a.get("chemical_formula_reduced"), a.get("_oqmd_spacegroup"), a.get("_oqmd_stability"))
    except Exception as e:
        print("OQMD", label, "ERR", type(e).__name__, str(e)[:100])
