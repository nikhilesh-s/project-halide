"""G4 convergence check for the proposed NaSbCl4 model: SCF at relaxed 60 Ry geometries with 80/640 Ry, k-spacing 0.17.
Run: discovery/.venv/bin/python -I overnight/make_converge_nasbcl4.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inputs as mi
mi.OUT = mi.ROOT / "overnight/night2/jobs"
for name, out in [("10c1_NaSbCl4_TP01", "05_NaSbCl4_TP01_CsAlCl4type_vcrelax.out"), ("10c2_NaSbCl4_agm005743441", "06_NaSbCl4_agm005743441_vcrelax.out"),
                  ("10c3_SbCl3_agm002205242", "07_SbCl3_agm002205242_vcrelax.out"), ("10c4_NaCl", "09_NaCl_agm003157609_vcrelax.out"),
                  ("10c5_Na3SbCl6", "10b_Na3SbCl6_agm017656109_vcrelax.out")]:
    s = mi.final_structure(mi.ROOT / "overnight/night2/jobs" / out)
    mi.write(name + "_80Ry_k017_scf", s, "scf", k=mi.kmesh(s, 0.17), ecut=80.0)
