"""Write night-3 QE inputs (queue in overnight/night3/jobs, overflow in overnight/night3/overflow). make_inputs.py is NOT edited:
its module source is exec'd with Br added to the hard-coded species tuple. Run from repo root."""
import json, types
from pathlib import Path
from pymatgen.core import Structure
ROOT = Path.cwd(); src = (ROOT/"overnight/make_inputs.py").read_text().replace('("Li", "Na", "Cl", "Bi", "Sb")', '("Li", "Na", "Cl", "Br", "Bi", "Sb")')
mi = types.ModuleType("mi"); mi.__file__ = str(ROOT/"overnight/make_inputs.py"); exec(compile(src, mi.__file__, "exec"), mi.__dict__)
mi.MASS["Br"] = 79.904; mi.UPF["Br"] = "Br.pbe-n-kjpaw_psl.1.0.0.UPF"
def cif(n):
    for d in ("gate3", "skeptic3"):
        p = ROOT/f"discovery/results/{d}/relax/{n}_d3.cif"
        if p.exists(): return Structure.from_file(p)
    raise FileNotFoundError(n)
# (file name, structure id, ecut, est. minutes)  -- estimates from night2/summary.json (12 at 60 Ry 31-73 min; LiBiCl4 12 at 80 Ry 71 min)
Q = [("01_Na3SbCl6_agm017656109", "agm017656109", 60, 40),
     ("02_LiSbBr4_TP01_hit", "TP-LiSbBr4-01", 80, 90), ("03_LiSbBr4_agm005739369_best_existing", "agm005739369", 80, 30),
     ("04_LiBr_rs_agm003157565_ref", "agm003157565", 80, 5), ("05_SbBr3_agm002139075_ref", "agm002139075", 80, 35),
     ("06_LiBiBr4_TP01_hit", "TP-LiBiBr4-01", 80, 90), ("07_LiBiBr4_agm005739442_best_existing", "agm005739442", 80, 30),
     ("08_BiBr3_agm002138960_ref", "agm002138960", 80, 35)]
O = [("09_NaSbBr4_TP01_hit", "TP-NaSbBr4-01", 80, 90), ("10_NaSbBr4_agm003733885_best_existing", "agm003733885", 80, 30),
     ("11_NaBr_agm003157587_ref", "agm003157587", 80, 5), ("12_Na3SbBr6_agm016963778", "agm016963778", 80, 60),
     ("13_LiSbCl4_TP03_hit", "TP-LiSbCl4-03", 80, 70), ("14_LiSbCl4_agm005453208_best_existing", "agm005453208", 80, 25),
     ("15_SbCl3_agm002142958_80Ry_ref", "agm002142958", 80, 20)]
meta = {}
for sub, lst in (("jobs", Q), ("overflow", O)):
    mi.OUT = ROOT/"overnight/night3"/sub; mi.OUT.mkdir(parents=True, exist_ok=True)
    for name, sid, ec, est in lst:
        s = cif(sid); mi.write(name + "_vcrelax", s, "vc-relax", ecut=float(ec)); meta[name] = dict(struct=sid, ecut=ec, est_min=est, queue=sub, natoms=len(s.get_primitive_structure()))
json.dump(meta, open(ROOT/"overnight/night3/queue.json", "w"), indent=1)
print("queue est h", sum(m["est_min"] for m in meta.values() if m["queue"] == "jobs")/60, "overflow est h", sum(m["est_min"] for m in meta.values() if m["queue"] == "overflow")/60)
