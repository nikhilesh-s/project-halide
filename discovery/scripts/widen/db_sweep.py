"""Widen Part A database sweep: every A-M-X chemical system in Alexandria PBE (OPTIMADE), Materials Project summary and
Materials Cloud MC3D (experimentally sourced structures). Resumable: one cached json per source x system.
Run: discovery/.venv/bin/python -I discovery/scripts/widen/db_sweep.py"""
import json, time, urllib.request, urllib.parse, itertools
from pathlib import Path
C = Path("discovery/results/widen/cache"); C.mkdir(parents=True, exist_ok=True)
A = ["Li", "Na", "K", "Rb", "Cs"]; M = ["Bi", "Sb", "In", "Ga", "Al", "Sc", "Y", "La", "Zn", "Mg", "Sn", "Pb"]; X = ["Cl", "Br", "I"]
KEY = Path.home().joinpath(".mp_api_key").read_text().strip()

def get(url, hdr=None):
    for k in range(4):
        try: return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=hdr or {}), timeout=60))
        except Exception as e: err = e; time.sleep(3*(k+1))
    raise err

def optimade(base, els, fields):
    url = base + "?" + urllib.parse.urlencode({"filter": "elements HAS ALL " + ",".join(f'"{e}"' for e in els) + " AND nelements=3",
                                              "page_limit": 100, "response_fields": fields})
    out = []
    while url:
        r = get(url); out += [{"id": d["id"], **d["attributes"]} for d in r["data"]]
        nx = (r.get("links") or {}).get("next"); url = nx.get("href") if isinstance(nx, dict) else nx
        time.sleep(0.5)
    return out

SRC = {
 "alex": lambda e: optimade("https://alexandria.icams.rub.de/pbe/v1/structures", e, "chemical_formula_reduced,nsites,_alexandria_hull_distance,_alexandria_space_group"),
 "mc3d": lambda e: optimade("https://optimade.materialscloud.org/main/mc3d-pbe-v1/structures", e, "chemical_formula_reduced,nsites,_mcloud_source_db,_mcloud_source_db_id"),
 "mp": lambda e: get("https://api.materialsproject.org/materials/summary/?" + urllib.parse.urlencode({"chemsys": "-".join(sorted(e)), "_limit": 200,
        "_fields": "material_id,formula_pretty,energy_above_hull,theoretical,symmetry,nsites,database_IDs"}), {"X-API-KEY": KEY, "User-Agent": "curl/8"})["data"],
}
if __name__ == "__main__":
    for m, a, x in itertools.product(M, A, X):  # M-major order: Bi, Sb systems first
        for s, fn in SRC.items():
            f = C/f"{s}_{a}-{m}-{x}.json"
            if f.exists(): continue
            try: json.dump(fn([a, m, x]), open(f, "w"))
            except Exception as e: print("FAIL", s, a, m, x, repr(e)[:120], flush=True)
        print(a, m, x, flush=True)
    print("DONE")
