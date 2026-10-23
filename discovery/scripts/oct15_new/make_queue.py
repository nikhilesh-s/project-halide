"""Write the 2026-10-15 new-candidate DFT inputs (front of overnight/m5_oct16/jobs; NaBiI4 long shot in m5_oct16/overflow,
not auto-run). Same approach as gate4/make_night4.py: make_inputs.py is NOT edited, its source is exec'd with I added.
Settings: PBE-D3, PSlibrary PAW, 80/640 Ry (bromide / iodide rule), k-spacing 0.22, vc-relax. Starting cells: rattled
MACE+D3 cells in discovery/results/oct15_new/relax. Run from repo root. TP-* = PROPOSED MODEL; nothing here is a claim."""
import json, shutil, types
from pathlib import Path
from pymatgen.core import Structure
ROOT = Path.cwd()
src = (ROOT/"overnight/make_inputs.py").read_text(); old = '("Li", "Na", "Cl", "Br", "In", "Sc", "Bi", "Sb")'; assert old in src
mi = types.ModuleType("mi"); mi.__file__ = str(ROOT/"overnight/make_inputs.py")
exec(compile(src.replace(old, '("Li", "Na", "Cl", "Br", "I", "In", "Sc", "Bi", "Sb")'), mi.__file__, "exec"), mi.__dict__)
mi.MASS["I"] = 126.90447; mi.UPF["I"] = "I.pbe-n-kjpaw_psl.1.0.0.UPF"   # same file and mass as night4 NaGaI4 jobs
cif = lambda n: Structure.from_file(ROOT/f"discovery/results/oct15_new/relax/{n}_d3.cif")
# (file name, structure id, est. M5 minutes). Done elsewhere, not repeated: NaBr agm003157587 80 Ry (air/15), BiBr3 R-3
# agm002138960 80 Ry (air/14), BiBr3 C2/m (m5_oct15/10i, tonight), NaI agm003157588 80 Ry (later/08, Air).
Q = [("00a_NaBiBr4_TP01_CsAlCl4type_hit", "TP-NaBiBr4-01", 90),
     ("00b_NaBiBr4_agm069058740_best_existing", "agm069058740", 30),
     ("00c_NaBiBr4_agm005739654_listed2", "agm005739654", 30),          # +3.1 above the best listed in MACE+D3: inside ML error
     ("00d_Na3BiBr6_agm053446239_ternary", "agm053446239", 50),         # Alexandria hull 0 (PBE): the competing ternary
     ("00e_NaSbI4_TP01_CsAlCl4type_hit", "TP-NaSbI4-01", 100),
     ("00f_NaSbI4_agm005743434_best_existing", "agm005743434", 35),
     ("00g_SbI3_agm003264389_ref", "agm003264389", 40)]                 # lowest SbI3 in MACE+D3 (Alexandria hull 0)
O = [("00h_NaBiI4_TP01_TlAlCl4type_hit", "TP-NaBiI4-01", 100),
     ("00i_NaBiI4_agm003740643_best_existing", "agm003740643", 35),
     ("00j_BiI3_agm003263390_ref", "agm003263390", 40)]                 # lowest BiI3 in MACE+D3
meta = {}
for sub, lst in (("jobs", Q), ("overflow", O)):
    mi.OUT = ROOT/"overnight/m5_oct16"/sub; mi.OUT.mkdir(parents=True, exist_ok=True)
    for name, sid, est in lst:
        assert not (mi.OUT/f"{name}_vcrelax.in").exists(), name
        s = cif(sid); mi.write(name + "_vcrelax", s, "vc-relax", ecut=80.0)
        meta[name] = dict(struct=sid, ecut=80, est_min=est, queue=sub, natoms=len(s.get_primitive_structure()), formula=s.composition.reduced_formula)
json.dump(meta, open(ROOT/"discovery/results/oct15_new/queue.json", "w"), indent=1)
S = {r["hit"]: r for r in json.load(open(ROOT/"discovery/results/oct15_new/summary.json"))}
for h in ("TP-NaBiBr4-01", "TP-NaSbI4-01", "TP-NaBiI4-01"):
    d = ROOT/"discovery/candidates"/h; assert d.exists()
    shutil.copy(ROOT/f"discovery/results/oct15_new/relax/{h}_d3.cif", d/"relaxed_mace_d3_rattled_oct15.cif")
    json.dump(dict(status="PROPOSED MODEL; MACE-MPA-0+D3 rattle check only (no DFT). DFT queued 2026-10-15, see discovery/results/new_candidates_oct15.md",
                   **{k: v for k, v in S[h].items() if k != "competitors"}), open(d/"oct15_check.json", "w"), indent=1, default=float)
print("jobs est h", sum(m["est_min"] for m in meta.values() if m["queue"] == "jobs")/60, "overflow est h", sum(m["est_min"] for m in meta.values() if m["queue"] == "overflow")/60)
