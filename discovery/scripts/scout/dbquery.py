"""Scout: count entries per chemical system in Alexandria PBE, OQMD (OPTIMADE) and Materials Project.
Usage: uv run --no-project python -I dbquery.py Na-Bi-Cl-Br Na-Bi-Cl ...
"""
import json, os, sys, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()

def get(url, headers=None, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "scout/0.1", **(headers or {})})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)

def optimade(base, els, fields, timeout):
    flt = "elements HAS ALL " + ",".join(f'"{e}"' for e in els) + f" AND nelements={len(els)}"
    url = base + "?" + urllib.parse.urlencode({"filter": flt, "response_fields": fields, "page_limit": 50})
    return flt, get(url, timeout=timeout)

def alex(els):
    try:
        flt, d = optimade("https://alexandria.icams.rub.de/pbe/v1/structures", els,
                          "chemical_formula_reduced,nsites,_alexandria_hull_distance,_alexandria_space_group", 90)
        rows = sorted((round(9 if x["attributes"].get("_alexandria_hull_distance") is None else x["attributes"]["_alexandria_hull_distance"], 4), x["attributes"]["chemical_formula_reduced"],
                 x["attributes"].get("nsites"), x["attributes"].get("_alexandria_space_group"), x["id"]) for x in d["data"])
        return f"ALEX n={d['meta'].get('data_returned', len(rows))} [{flt}] " + "; ".join(map(str, rows[:8]))
    except Exception as e:
        return f"ALEX ERROR {type(e).__name__}: {str(e)[:100]}"

def oqmd(els):
    try:
        flt, d = optimade("https://oqmd.org/optimade/structures", els,
                          "chemical_formula_reduced,nsites,_oqmd_stability,_oqmd_entry_id,_oqmd_spacegroup", 420)
        rows = sorted((round(9 if x["attributes"].get("_oqmd_stability") is None else x["attributes"]["_oqmd_stability"], 4), x["attributes"]["chemical_formula_reduced"],
                 x["attributes"].get("nsites"), x["attributes"].get("_oqmd_spacegroup"), x["attributes"].get("_oqmd_entry_id")) for x in d["data"])
        return f"OQMD n={d['meta'].get('data_returned', len(rows))} [{flt}] " + "; ".join(map(str, rows[:8]))
    except Exception as e:
        return f"OQMD ERROR {type(e).__name__}: {str(e)[:100]}"

def mp(els):
    try:
        cs = "-".join(sorted(els))
        d = get("https://api.materialsproject.org/materials/summary/?" + urllib.parse.urlencode(
            {"chemsys": cs, "_fields": "material_id,formula_pretty,energy_above_hull,nsites,theoretical,symmetry", "_limit": 50}),
            {"X-API-KEY": KEY}, 60)
        rows = sorted((round(x.get("energy_above_hull") or 0, 4), x["formula_pretty"], x.get("nsites"),
                 (x.get("symmetry") or {}).get("symbol"), "theo" if x.get("theoretical") else "EXPT", x["material_id"]) for x in d["data"])
        return f"MP   n={len(rows)} [chemsys={cs}] " + "; ".join(map(str, rows[:8]))
    except Exception as e:
        return f"MP   ERROR {type(e).__name__}: {str(e)[:100]}"

if __name__ == "__main__":
    dbs = [alex, mp] + ([] if os.environ.get("NO_OQMD") else [oqmd])
    systems = [s.split("-") for s in sys.argv[1:]]
    with ThreadPoolExecutor(12) as ex:
        futs = [(s, [ex.submit(f, s) for f in dbs]) for s in systems]
        for s, fs in futs:
            print("==", "-".join(s))
            for f in fs:
                print("  ", f.result())
