"""Write night-4 QE inputs (overnight/night4/jobs, overflow in overnight/night4/overflow). make_inputs.py is NOT edited:
its source is exec'd with the species tuple extended, and MASS/UPF patched for the new elements. Run from repo root.
Cutoff rule: ecutwfc = max(60, highest PSlibrary 'suggested minimum cutoff' among the elements, rounded up to 60/80/100);
Li (103 Ry suggested) is capped at the top tier, 100 Ry. One cutoff per reaction."""
import json, shutil, types
from pathlib import Path
from pymatgen.core import Structure
ROOT = Path.cwd()
SP = '("Li", "Na", "Cl", "Br", "I", "Al", "Ga", "In", "Sc", "Y", "Bi", "Sb")'
src = (ROOT/"overnight/make_inputs.py").read_text(); assert '("Li", "Na", "Cl", "Bi", "Sb")' in src
mi = types.ModuleType("mi"); mi.__file__ = str(ROOT/"overnight/make_inputs.py")
exec(compile(src.replace('("Li", "Na", "Cl", "Bi", "Sb")', SP), mi.__file__, "exec"), mi.__dict__)
mi.MASS.update(Br=79.904, I=126.90447, Al=26.98154, Ga=69.723, In=114.818, Sc=44.95591, Y=88.90584)
mi.UPF.update(Br="Br.pbe-n-kjpaw_psl.1.0.0.UPF", I="I.pbe-n-kjpaw_psl.1.0.0.UPF", Al="Al.pbe-n-kjpaw_psl.1.0.0.UPF", Ga="Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
              In="In.pbe-dn-kjpaw_psl.1.0.0.UPF", Sc="Sc.pbe-spn-kjpaw_psl.1.0.0.UPF", Y="Y.pbe-spn-kjpaw_psl.1.0.0.UPF")
SUG = dict(Li=103, Na=66, Cl=45, Br=49, I=36, Al=29, Ga=60, In=41, Sc=50, Y=38, Bi=45, Sb=27)   # Ry, read from the UPF headers
def ecut(chem_els): return next(t for t in (60, 80, 100) if t >= min(max(SUG[e] for e in chem_els), 100))
assert ecut(["Li", "In", "Cl"]) == 100 and ecut(["Na", "Ga", "I"]) == 80 and ecut(["Cl", "Bi"]) == 60
def cif(n):
    for d in ("gate4", "gate3", "skeptic3"):
        p = ROOT/f"discovery/results/{d}/relax/{n}_d3.cif"
        if p.exists(): return Structure.from_file(p)
    raise FileNotFoundError(n)
# (file name, structure id, reaction elements, est. minutes). Estimates: night2 12 atoms at 60 Ry 31-73 min; x2 at 80 Ry, x3.5 at 100 Ry.
LIC, LISB, LIYB, NAGI = ["Li", "In", "Cl"], ["Li", "Sc", "Br"], ["Li", "Y", "Br"], ["Na", "Ga", "I"]
Q = [("01_LiInCl4_TP01_hit", "TP-LiInCl4-01", LIC, 150), ("02_LiInCl4_agm003740637_best_existing", "agm003740637", LIC, 50),
     ("03_LiCl_wz_agm002157414_100Ry_ref", "agm002157414", LIC, 15), ("04_InCl3_agm002142774_ref", "agm002142774", LIC, 70),
     ("05_LiScBr4_TP01_hit", "TP-LiScBr4-01", LISB, 150), ("06_LiScBr4_agm005739333_best_existing", "agm005739333", LISB, 50),
     ("07_LiBr_agm002091979_100Ry_ref", "agm002091979", LISB, 15), ("08_ScBr3_agm002184240_ref", "agm002184240", LISB, 70),
     ("09_LiSc3Br8_agm006022042_ternary", "agm006022042", LISB, 150),
     ("10_LiYBr4_TP02_hit", "TP-LiYBr4-02", LIYB, 150), ("11_LiYBr4_agm003742437_best_existing", "agm003742437", LIYB, 50),
     ("12_YBr3_agm003297005_ref", "agm003297005", LIYB, 70),
     ("13_NaGaI4_TP02_hit", "TP-NaGaI4-02", NAGI, 90), ("14_NaGaI4_agm003736187_best_existing", "agm003736187", NAGI, 90),
     ("15_NaGaI4_agm005755304_alex_hull0", "agm005755304", NAGI, 30), ("16_NaI_agm003157588_ref", "agm003157588", NAGI, 5),
     ("17_GaI3_agm003220383_ref", "agm003220383", NAGI, 150)]
O = [("18_NaGaI3_agm005838860_ternary", "agm005838860", NAGI, 200)]
meta = {}
for sub, lst in (("jobs", Q), ("overflow", O)):
    mi.OUT = ROOT/"overnight/night4"/sub; mi.OUT.mkdir(parents=True, exist_ok=True)
    for name, sid, els, est in lst:
        s = cif(sid); ec = ecut(els); mi.write(name + "_vcrelax", s, "vc-relax", ecut=float(ec))
        meta[name] = dict(struct=sid, ecut=ec, est_min=est, queue=sub, natoms=len(s.get_primitive_structure()), formula=s.composition.reduced_formula)
json.dump(meta, open(ROOT/"overnight/night4/queue.json", "w"), indent=1)
S = {r["hit"]: r for r in json.load(open(ROOT/"discovery/results/gate4/summary.json"))}
for h in ("TP-LiInCl4-01", "TP-LiScBr4-01", "TP-LiYBr4-02", "TP-NaGaI4-02"):   # kept hits: PROPOSED MODELS
    d = ROOT/"discovery/candidates"/h; assert d.exists()
    shutil.copy(ROOT/f"discovery/results/gate4/relax/{h}_d3.cif", d/"relaxed_mace_d3_rattled_gate4.cif")
    json.dump(dict(status="PROPOSED MODEL; passed gate 4 (MACE-MPA-0+D3 rattle check only, no DFT)", **{k: v for k, v in S[h].items() if k != "competitors"}),
              open(d/"gate4_check.json", "w"), indent=1, default=float)
print("queue est h", sum(m["est_min"] for m in meta.values() if m["queue"] == "jobs")/60, "overflow est h", sum(m["est_min"] for m in meta.values() if m["queue"] == "overflow")/60)
