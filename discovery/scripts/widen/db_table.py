"""Summarise the db_sweep cache: per composition, Alexandria (n entries, min hull, min nsites, sg), MP (ehull, theoretical, ICSD ids), MC3D source.
Run: discovery/.venv/bin/python -I discovery/scripts/widen/db_table.py [M ...]   -> discovery/results/widen/db_table.csv"""
import json, sys, csv, glob, re
from pathlib import Path
from pymatgen.core import Composition
C = Path("discovery/results/widen/cache"); rows = {}
red = lambda f: Composition(f).reduced_formula
for f in sorted(C.glob("*.json")):
    src, cs = f.stem.split("_"); d = json.load(open(f))
    for e in d:
        try:
            if src == "alex":
                r = rows.setdefault((cs, red(e["chemical_formula_reduced"])), {}); h = e["_alexandria_hull_distance"]
                r["alex_n"] = r.get("alex_n", 0) + 1
                if h < r.get("alex_hull", 9): r.update(alex_hull=h, alex_nsites_best=e["nsites"], alex_sg_best=e["_alexandria_space_group"], alex_id=e["id"])
                r["alex_min_nsites"] = min(r.get("alex_min_nsites", 999), e["nsites"])
            elif src == "mp":
                r = rows.setdefault((cs, red(e["formula_pretty"])), {}); ic = (e.get("database_IDs") or {}).get("icsd", [])
                r["mp"] = r.get("mp", "") + f'{e["material_id"]}:{e["symmetry"]["symbol"]}:n{e["nsites"]}:eh{1000*e["energy_above_hull"]:.0f}:theo={e["theoretical"]}:icsd{len(ic)}; '
            else:
                r = rows.setdefault((cs, red(e["chemical_formula_reduced"])), {}); r["mc3d"] = r.get("mc3d", "") + f'{e.get("_mcloud_source_db")}:{e.get("_mcloud_source_db_id")}:n{e["nsites"]}; '
        except Exception as ex: print("skip", f.name, repr(ex)[:80])
cols = ["chemsys", "formula", "alex_n", "alex_hull", "alex_nsites_best", "alex_min_nsites", "alex_sg_best", "alex_id", "mp", "mc3d"]
with open("discovery/results/widen/db_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, cols); w.writeheader()
    for (cs, fm), r in sorted(rows.items()): w.writerow(dict(chemsys=cs, formula=fm, **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
want = sys.argv[1:]
for (cs, fm), r in sorted(rows.items()):
    if want and cs.split("-")[1] not in want: continue
    exp = ("EXP(mc3d)" if r.get("mc3d") else "") + (" MPexp" if "theo=False" in r.get("mp", "") else "") + (" MPtheo" if "theo=True" in r.get("mp", "") else "")
    if r.get("alex_hull", 9) > 0.05 and not exp.strip(): continue
    print(f'{cs:9s} {fm:10s} alex n={r.get("alex_n",0):2d} hull={1000*r.get("alex_hull",9):6.1f} nbest={r.get("alex_nsites_best","-")} nmin={r.get("alex_min_nsites","-")} sg={r.get("alex_sg_best","-")} | {exp}')
