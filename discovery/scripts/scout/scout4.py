"""Cycle-4 scout: MP summary + COD element search for ternary chemsys given as args (e.g. Li-Al-Br)."""
import json, os, sys, urllib.request, urllib.parse

key = open(os.path.expanduser("~/.mp_api_key")).read().strip()
F = "material_id,formula_pretty,energy_above_hull,theoretical,symmetry,nsites,volume,database_IDs"


def get(u, h=None, t=60):
    return json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "curl/8", **(h or {})}), timeout=t))


for cs in sys.argv[1:]:
    els = cs.split("-")
    try:
        d = get("https://api.materialsproject.org/materials/summary/?" + urllib.parse.urlencode(
            {"chemsys": "-".join(sorted(els)), "_fields": F, "_limit": 200}), {"X-API-KEY": key, "accept": "application/json"})["data"]
        print("MP", cs, len(d))
        for x in sorted(d, key=lambda x: x["formula_pretty"]):
            print("  ", x["material_id"], x["formula_pretty"], "Ehull=%.3f" % x["energy_above_hull"], "theor=%s" % x["theoretical"],
                  x["symmetry"]["symbol"], "nsites=%d" % x["nsites"], "V=%.1f" % x["volume"], (x.get("database_IDs") or {}).get("icsd", [])[:6])
    except Exception as e:
        print("MP", cs, "ERR", type(e).__name__, str(e)[:100])
    try:
        q = {"el%d" % (i + 1): e for i, e in enumerate(els)}
        q.update(strictmin=len(els), strictmax=len(els), format="json")
        d = get("https://www.crystallography.net/cod/result?" + urllib.parse.urlencode(q), t=25)
        print("COD", cs, len(d))
        for x in d:
            print("  ", x.get("file"), x.get("formula"), x.get("sg"), "Z=%s" % x.get("Z"), x.get("a"), x.get("b"), x.get("c"), x.get("beta"),
                  "|", (x.get("authors") or "")[:40], x.get("journal"), x.get("year"), x.get("doi"))
    except Exception as e:
        print("COD", cs, "ERR", type(e).__name__, str(e)[:100])
