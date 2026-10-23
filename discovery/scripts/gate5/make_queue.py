"""Gate 5 step 2: rewrite the remaining DFT queue at ONE cutoff per candidate set. Run from repo root:
  discovery/.venv/bin/python -I discovery/scripts/gate5/make_queue.py
overnight/night4/jobs (about 7 h), overnight/night5/jobs (about 7 h), overnight/later/jobs (not scheduled; `run_night.py later`).
The old 100 Ry night-4 inputs are moved to overnight/later/superseded_100Ry/, not deleted. make_inputs.py is NOT edited
(same exec/monkey-patch as gate4/make_night4.py).
Cutoff: 80/640 Ry for every Li or Na+Ga job, i.e. the cutoff already used for all Li work (night2, night3, air_day), so
LiCl (night2 12-13) and LiBr rs (night3 04) are reused. Li.pbe-s-kjpaw suggests 103 Ry; a 100 Ry + k 0.17 SCF check on the
relaxed cells is the G4 step for any candidate that passes G2-G3.
Estimates (min, M5): night2 timings at 60 Ry (12 atoms 31-73, 10 atoms 37, 8 atoms 10-28, 6 atoms 8-21) x2 for 80 Ry."""
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
ECUT = 80.0
def cif(n):
    for d in ("gate5", "gate4", "gate3", "skeptic3"):
        p = ROOT/f"discovery/results/{d}/relax/{n}_d3.cif"
        if p.exists(): return Structure.from_file(p)
    raise FileNotFoundError(n)
# (name without number, structure id, est. minutes, judges candidate)
PLAN = {
 "night4": [("Li3BiCl6_agm005963711_ternary", "agm005963711", 75, "TP-LiBiCl4-01"), ("Li3SbBr6_agm006027634_ternary", "agm006027634", 75, "TP-LiSbBr4-01"),
            ("Li3BiBr6_agm005961056_ternary", "agm005961056", 75, "TP-LiBiBr4-01"), ("Li3SbCl6_agm006027845_ternary", "agm006027845", 75, "TP-LiSbCl4-03"),
            ("LiInCl4_TP01_hit", "TP-LiInCl4-01", 100, "TP-LiInCl4-01"), ("LiInCl4_agm003740637_best_existing", "agm003740637", 30, "TP-LiInCl4-01")],
 "night5": [("InCl3_agm002142774_ref", "agm002142774", 45, "TP-LiInCl4-01"), ("Li3InCl6_agm003569151_ternary", "agm003569151", 75, "TP-LiInCl4-01"),
            ("LiScBr4_TP01_hit", "TP-LiScBr4-01", 100, "TP-LiScBr4-01"), ("LiScBr4_agm005739333_best_existing", "agm005739333", 30, "TP-LiScBr4-01"),
            ("LiBr_agm002091979_ref", "agm002091979", 10, "TP-LiScBr4-01, TP-LiYBr4-02"), ("ScBr3_agm002184240_ref", "agm002184240", 45, "TP-LiScBr4-01"),
            ("LiSc3Br8_agm006022042_ternary", "agm006022042", 100, "TP-LiScBr4-01")],
 "later":  [("LiYBr4_TP02_hit", "TP-LiYBr4-02", 100, "TP-LiYBr4-02"), ("LiYBr4_agm003742437_best_existing", "agm003742437", 30, "TP-LiYBr4-02"),
            ("YBr3_agm003297005_ref", "agm003297005", 45, "TP-LiYBr4-02"), ("Li3YBr6_agm003224330_ternary", "agm003224330", 300, "TP-LiYBr4-02"),
            ("NaGaI4_TP02_hit", "TP-NaGaI4-02", 90, "TP-NaGaI4-02"), ("NaGaI4_agm003736187_best_existing", "agm003736187", 90, "TP-NaGaI4-02"),
            ("NaGaI4_agm005755304_alex_hull0", "agm005755304", 30, "TP-NaGaI4-02"), ("NaI_agm003157588_ref", "agm003157588", 5, "TP-NaGaI4-02"),
            ("GaI3_agm003220383_ref", "agm003220383", 150, "TP-NaGaI4-02"), ("NaGaI3_agm005838860_ternary", "agm005838860", 200, "TP-NaGaI4-02 (Ga(II): off the NaI-GaI3 line)"),
            ("Li3ScBr6_agm003568530_ternary_optional", "agm003568530", 300, "TP-LiScBr4-01 (optional; +67 meV in Alexandria)")]}
old = ROOT/"overnight/later/superseded_100Ry"; old.mkdir(parents=True, exist_ok=True)
for q in PLAN:   # clear previous inputs (never outputs): night-4 originals are kept under later/superseded_100Ry
    for f in list((ROOT/"overnight"/q).glob("jobs/*.in")) + list((ROOT/"overnight"/q).glob("overflow/*.in")):
        assert not f.with_suffix(".out").exists(), f"refusing to replace an input that has output: {f}"
        if "ecutwfc=100" in f.read_text() or f.parent.name == "overflow": shutil.move(str(f), str(old/f.name))
        else: f.unlink()   # regenerated below by this script
if (ROOT/"overnight/night4/overflow").exists() and not any((ROOT/"overnight/night4/overflow").iterdir()): (ROOT/"overnight/night4/overflow").rmdir()
if (ROOT/"overnight/night4/queue.json").exists() and not (old/"queue.json").exists(): shutil.copy(ROOT/"overnight/night4/queue.json", old/"queue.json")
allmeta = {}
for q, lst in PLAN.items():
    mi.OUT = ROOT/"overnight"/q/"jobs"; mi.OUT.mkdir(parents=True, exist_ok=True); meta = {}
    for n, (name, sid, est, cand) in enumerate(lst, 1):
        s = cif(sid); full = f"{n:02d}_{name}_vcrelax"; mi.write(full, s, "vc-relax", ecut=ECUT)
        t = (mi.OUT/f"{full}.in").read_text(); assert f"prefix='{full}'" in t and f"outdir='__TMP__/{full}'" in t and "ecutwfc=80.0" in t
        meta[full] = dict(struct=sid, formula=s.composition.reduced_formula, natoms=len(s.get_primitive_structure()), ecut=ECUT, est_min=est, judges=cand, queue=q)
    json.dump(meta, open(ROOT/"overnight"/q/"queue.json", "w"), indent=1); allmeta.update(meta)
    print(q, len(lst), "jobs, est h %.1f" % (sum(x[2] for x in lst)/60))
json.dump(allmeta, open(ROOT/"discovery/results/gate5/queue.json", "w"), indent=1)
