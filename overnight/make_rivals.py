"""DFT inputs for the two NaBiCl4 scan models not yet tested by DFT (full-scan names TP-NaBiCl4-01 KGaI4-type
and TP-NaBiCl4-03 NaGaF4-type). NOTE: the DFT-tested model called "TP-NaBiCl4-01" in EVIDENCE/SCOREBOARD is the
full scan's TP-NaBiCl4-02 (KBiF4-type, agm077914577); the full scan renumbered candidates.
Run: discovery/.venv/bin/python -I overnight/make_rivals.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inputs as mi
from pymatgen.core import Structure
mi.OUT = mi.ROOT / "overnight/night2/jobs"
for name, cand in [("10a1_NaBiCl4_TPscan01_KGaI4type_RIVAL", "TP-NaBiCl4-01"), ("10a2_NaBiCl4_TPscan03_NaGaF4type_RIVAL", "TP-NaBiCl4-03")]:
    mi.write(name + "_vcrelax", Structure.from_file(mi.ROOT / f"discovery/candidates/{cand}/relaxed_mace_d3.cif"), "vc-relax", ecut=60.0)
