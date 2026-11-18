"""Tabulate skeptic results from discovery/results/skeptic/relax/*.json. Run from repo root."""
import json, glob, math
R = {json.load(open(f))["name"]: json.load(open(f)) for f in glob.glob("discovery/results/skeptic/relax/*.json")}
A = json.load(open("discovery/results/skeptic/alex.json")); AE = {e["id"]: e for v in A.values() for e in v}
pbe = lambda i: AE[i]["_alexandria_energy"]/AE[i]["nsites"]
nc = [k for k, r in R.items() if r["relaxed"] and not r["conv"]]; print("NOT CONVERGED:", nc)
for form, ids in [("NaBiBr4", [e["id"] for e in A["BiBr4Na"]]), ("NaBiCl4", [e["id"] for e in A["BiCl4Na"]])]:
    ref = min(ids, key=pbe); print(f"\n{form}: energies relative to {ref} (meV/atom); V in A^3/atom")
    print("id | PBE(Alex) | hull | MACE@PBEgeom | MACE relaxed | MACE+D3 relaxed | V_PBE | V_MACE | V_MACE+D3")
    rows = ids + (["model"] if form == "NaBiBr4" else [])
    for i in rows:
        g = lambda t: 1000*(R[f"{i}_{t}"]["E"]-R[f"{ref}_{t}"]["E"])
        if i == "model": print(f"model | - | - | - | {g('nod'):+.1f} | {g('d3'):+.1f} | - | {R['model_nod']['V']:.2f} | {R['model_d3']['V']:.2f}")
        else: print(f"{i} | {1000*(pbe(i)-pbe(ref)):+.1f} | {1000*AE[i]['_alexandria_hull_distance']:.1f} | {g('sp'):+.1f} | {g('nod'):+.1f} | {g('d3'):+.1f} | {R[i+'_sp']['V']:.2f} | {R[i+'_nod']['V']:.2f} | {R[i+'_d3']['V']:.2f}")
print("\nAll NaBiBr4 candidates incl. prototype decorations, rel. to agm069058740 (meV/atom)")
for t in ("nod", "d3"):
    c = sorted(((1000*(r["E"]-R[f"agm069058740_{t}"]["E"]), k, r["V"], r["conv"]) for k, r in R.items() if k.endswith("_"+t) and r["formula"] == "NaBiBr4"))
    for e, k, v, cv in c: print(f"  {t} {e:+7.1f}  V={v:.2f} conv={cv} {k}")
print("\nFixed volume (no D3), E rel. to agm069058740 at same V (meV/atom); abs E of agm069058740")
names = ["model", "agm069058740", "agm005739654", "agm003733888"]
for v in (30, 33, 36, 39, 42):
    print(f"  V={v}: " + " ".join(f"{n}={1000*(R[f'fixV{v}_{n}']['E']-R[f'fixV{v}_agm069058740']['E']):+.1f}({'c' if R[f'fixV{v}_{n}']['conv'] else 'NC'})" for n in names), f"E0={R[f'fixV{v}_agm069058740']['E']:.4f}")
print("\nDecomposition (meV/atom of NaBiBr4; negative = compound stable)")
nabr, bibr3, n3 = "agm003157587", "agm002138960", "agm053446239"
def rx(E, x):  # x = E/atom of NaBiBr4
    return 1000*(x-(2*E(nabr)+4*E(bibr3))/6), 1000*(x-(10*E(n3)+8*E(bibr3))/18)
for i in [e["id"] for e in A["BiBr4Na"]]:
    print(f"  PBE Alexandria {i}: vs NaBr+BiBr3 {rx(pbe, pbe(i))[0]:+.1f}; vs Na3BiBr6+2BiBr3 {rx(pbe, pbe(i))[1]:+.1f} (Alex hull {1000*AE[i]['_alexandria_hull_distance']:.1f})")
for t in ("sp", "nod", "d3"):
    E = lambda i: R[f"{i}_{t}"]["E"]
    for i in [e["id"] for e in A["BiBr4Na"]] + ([] if t == "sp" else ["model"]):
        a, b = rx(E, E(i)); print(f"  MACE {t} {i}: vs NaBr+BiBr3 {a:+.1f}; vs Na3BiBr6+2BiBr3 {b:+.1f}")
    print(f"    Na3BiBr6 vs 3NaBr+BiBr3 ({t}): {1000*(E(n3)-(6*E(nabr)+4*E(bibr3))/10):+.1f} meV/atom" )
print(f"    Na3BiBr6 vs 3NaBr+BiBr3 (PBE): {1000*(pbe(n3)-(6*pbe(nabr)+4*pbe(bibr3))/10):+.1f}")
for t in ("nod", "d3"): print(t, "binary volumes", {i: (round(R[f'{i}_sp']['V'], 2), round(R[f'{i}_{t}']['V'], 2)) for i in (nabr, bibr3, n3)})
kT = 8.617333e-5*300
print("\nIdeal mixing entropy term -T*S_conf at 300 K, meV per atom (4 anion sites per 6 atoms):")
for x in (1, 2, 3):
    p = x/4; s = -(p*math.log(p)+(1-p)*math.log(1-p)); print(f"  x={x}: {-1000*kT*s*4/6:.1f} meV/atom ({-1000*kT*s:.1f} per anion site)")
