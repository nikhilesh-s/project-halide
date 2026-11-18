"""Attack 5: redo the 42 mixed Cl/Br orderings with MACE-MPA-0 + D3(BJ), starting from the cycle-2 relaxed cells."""
import sys, json, glob; sys.path.insert(0, "discovery/scripts/skeptic")
from pathlib import Path
from pymatgen.core import Structure
from run import job, OUT
R = {}
for d in sorted(glob.glob("discovery/results/mace/NBCB-*")):
    i = Path(d).name; R[i] = job(f"mix_{i}_d3", Structure.from_file(f"{d}/relaxed.cif"), True)
ld = lambda n: json.load(open(OUT/f"{n}.json"))["E"]
ECl, EBr = R["NBCB-0-parent"]["E"], ld("model_d3")            # same-prototype end members
ECl0, EBr0 = ld("agm088474086_d3"), ld("agm069058740_d3")      # lowest MACE+D3 end members found
rows = []
for i, r in R.items():
    x = int(i.split("-")[1])
    if 0 < x < 4: rows.append((1000*(r["E"]-((4-x)*ECl+x*EBr)/4), 1000*(r["E"]-((4-x)*ECl0+x*EBr0)/4), x, i, r["conv"]))
rows.sort()
print("parent-vs-agm005743614_d3 check (meV):", 1000*(ECl-ld("agm005743614_d3")))
for x in (1, 2, 3):
    s = [r for r in rows if r[2] == x]; print(f"x={x}: n={len(s)} E_mix same-prototype min {s[0][0]:+.2f} ({s[0][3]}) max {s[-1][0]:+.2f}; vs lowest end members min {min(r[1] for r in s):+.2f}")
print("all:", [(r[3], round(r[0], 2)) for r in rows[:6]], "nonconv:", [r[3] for r in rows if not r[4]])
