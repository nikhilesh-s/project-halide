"""Energies from overnight/night1+night2 (QE PBE-D3, last '!' = fresh-basis final SCF). Prints meV/atom tables."""
import json, re, glob
from pathlib import Path
from pymatgen.core import Structure, Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
RY = 13605.693123
def parse(out):
    t = Path(out).read_text(errors="replace")
    e = float(re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", t, re.M)[-1])
    blk = t[t.rindex("CELL_PARAMETERS"):].splitlines()
    m = re.search(r"alat=\s*([\d.]+)", blk[0]); sc = float(m.group(1)) * 0.529177210903 if m else 1.0
    cell = [[float(x) * sc for x in l.split()] for l in blk[1:4]]
    i = next(n for n, l in enumerate(blk) if l.startswith("ATOMIC_POSITIONS"))
    sp, xyz = [], []
    for l in blk[i + 1:]:
        w = l.split()
        if len(w) != 4 or not w[0].isalpha(): break
        sp.append(w[0]); xyz.append([float(x) for x in w[1:]])
    s = Structure(cell, sp, xyz)
    return dict(e=e, s=s, n=len(s), conv="bfgs converged" in t, f=s.composition.reduced_formula,
                fu=s.composition.get_reduced_composition_and_factor()[1], sg=SpacegroupAnalyzer(s, 0.05).get_space_group_symbol(), v=s.volume / len(s))
R = {Path(f).stem.replace("_vcrelax", ""): parse(f) for f in sorted(glob.glob("overnight/night[12]/jobs/*vcrelax.out"))}
for k, r in R.items(): print(f"{k:42s} {r['f']:9s} n={r['n']:2d} {r['sg']:9s} V/at={r['v']:.2f} conv={r['conv']} E/fu={r['e']/r['fu']:.6f}")
epf = lambda k: R[k]["e"] / R[k]["fu"]
best = lambda form: min((k for k in R if R[k]["f"] == form and R[k]["conv"]), key=epf)
print("\nNaBiCl4 polymorph ranking (meV/atom vs re-relaxed agm005743614):")
ref = epf("00_NaBiCl4_agm005743614_rerelax")
for k in sorted((k for k in R if R[k]["f"] == "NaBiCl4"), key=epf): print(f"  {k:42s} {(epf(k)-ref)*RY/6:8.2f}  {R[k]['sg']}")
nacl, bicl3, sbcl3, na3 = epf(best("NaCl")), epf(best("BiCl3")), epf(best("SbCl3")), epf(best("Na3BiCl6"))
print("references:", best("NaCl"), best("BiCl3"), best("SbCl3"))
for k in [k for k in R if R[k]["f"] == "NaBiCl4"]:
    print(f"  {k:42s} vs NaCl+BiCl3 {(epf(k)-nacl-bicl3)*RY/6:7.2f} | vs (Na3BiCl6+2BiCl3)/3 {(3*epf(k)-na3-2*bicl3)*RY/18:7.2f}")
for k in [k for k in R if R[k]["f"] == "NaSbCl4"]:
    print(f"  {k:42s} vs NaCl+SbCl3 {(epf(k)-nacl-sbcl3)*RY/6:7.2f}")
sm = StructureMatcher()
nb = [k for k in R if R[k]["f"] == "NaBiCl4"]
print("TP01 matches any other NaBiCl4:", [k for k in nb if "TP01" not in k and sm.fit(R["01_NaBiCl4_TP01_KBiF4type"]["s"], R[k]["s"])])
print("NaSbCl4 TP01 == agm005743441:", sm.fit(R["05_NaSbCl4_TP01_CsAlCl4type"]["s"], R["06_NaSbCl4_agm005743441"]["s"]))
