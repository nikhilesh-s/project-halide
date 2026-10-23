"""Sr3Ta2O5N2 novelty check: Crossref + Alexandria OPTIMADE + Materials Project. No DFT, no ML.
Usage (repo root): discovery/.venv/bin/python -I discovery/scripts/sr3ta2o5n2/novelty_query.py
"""
import json, os, time, urllib.parse, urllib.request

KEY = open(os.path.expanduser("~/.mp_api_key")).read().strip()   # never printed

def get(url, headers=None):
    req = urllib.request.Request(url, headers={"User-Agent": "halide-scout/0.1", **(headers or {})})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)

CROSSREF = ["Sr3Ta2O5N2", "Sr3Ta2N2O5", "strontium tantalum oxynitride Ruddlesden-Popper n=2",
            "Sr3Nb2O5N2", "Anion ordering Ruddlesden-Popper strontium niobium oxynitrides"]
FORMULAS = {"Sr3Ta2O5N2": "N2O5Sr3Ta2", "Sr3Nb2O5N2": "Nb2N2O5Sr3", "Sr2TaO3N": "NO3Sr2Ta", "SrTaO2N": "NO2SrTa",
            "Ta3N5": "N5Ta3", "TaON": "NOTa", "Sr4Ta2O9": "O9Sr4Ta2", "Sr5Ta4O15": "O15Sr5Ta4",
            "SrO": "OSr", "Sr3N2": "N2Sr3", "Sr2N": "NSr2"}

def crossref(q):
    d = get("https://api.crossref.org/works?" + urllib.parse.urlencode({"query.bibliographic": q, "rows": 20, "select": "title,DOI,issued"}))
    for it in d["message"]["items"]:
        t = " ".join((it.get("title") or ["-"])[0].split())
        if any(k in t for k in ("Sr", "tantal", "Tantal", "oxynitr", "Oxynitr", "Ruddlesden")):
            print(f"   {it['DOI']}  {t[:140]}")

def alex(f):
    d = get("https://alexandria.icams.rub.de/pbe/v1/structures?" + urllib.parse.urlencode(
        {"filter": f'chemical_formula_reduced="{f}"', "response_fields": "nsites,_alexandria_hull_distance,_alexandria_space_group"}))
    return [(x["id"], x["attributes"].get("nsites"), x["attributes"].get("_alexandria_space_group"),
             x["attributes"].get("_alexandria_hull_distance")) for x in d["data"]]

def mp(name):
    d = get("https://api.materialsproject.org/materials/summary/?" + urllib.parse.urlencode(
        {"formula": name, "_fields": "material_id,nsites,energy_above_hull,theoretical,symmetry", "_limit": 20}), {"X-API-KEY": KEY})
    return [(x["material_id"], x.get("nsites"), (x.get("symmetry") or {}).get("symbol"),
             round(x.get("energy_above_hull") or 0, 4), "theo" if x.get("theoretical") else "EXPT") for x in d["data"]]

if __name__ == "__main__":
    for q in CROSSREF:
        print("CROSSREF", q); time.sleep(3)
        try: crossref(q)
        except Exception as e: print("   ERROR", type(e).__name__, str(e)[:80])
    for name, f in FORMULAS.items():
        for tag, fn, arg in (("ALEX", alex, f), ("MP", mp, name)):
            try: print(f"{tag:4} {name:11}", fn(arg))
            except Exception as e: print(f"{tag:4} {name:11} ERROR {type(e).__name__}: {str(e)[:80]}")
