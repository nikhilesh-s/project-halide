"""Convergence check (Track C1): single-point SCF at the relaxed 60 Ry geometries, with 80/640 Ry and a denser
k-mesh (spacing 0.17 1/A), for the NaBiCl4 polymorphs and their competitors. Run from repo root:
discovery/.venv/bin/python -I overnight/make_converge.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inputs as mi
ROOT = mi.ROOT
mi.OUT = ROOT / "overnight/converge_air/jobs"; mi.OUT.mkdir(parents=True, exist_ok=True)
SRC = [("01_NaBiCl4_TP01", "night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.out"),
       ("02_NaBiCl4_agm088474086", "night1/jobs/01_NaBiCl4_agm088474086_vcrelax.out"),
       ("03_NaCl", "night2/jobs/09_NaCl_agm003157609_vcrelax.out"),
       ("04_BiCl3", "night2/jobs/10_BiCl3_agm002142654_vcrelax.out"),
       ("05_Na3BiCl6", "night2/jobs/02_Na3BiCl6_agm067313194_vcrelax.out"),
       ("06_NaBiCl4_agm084661913", "night1/jobs/02_NaBiCl4_agm084661913_vcrelax.out"),
       ("07_NaBiCl4_agm005743614", "night2/jobs/00_NaBiCl4_agm005743614_rerelax_vcrelax.out")]
for name, out in SRC:
    s = mi.final_structure(ROOT / "overnight" / out)
    mi.write(name + "_80Ry_k017_scf", s, "scf", k=mi.kmesh(s, 0.17), ecut=80.0)
