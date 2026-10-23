"""Inputs for the electrochemical stability window and band gap of the top two (TP-NaSbBr4-01, TP-LiSbBr4-01).
Elemental references Na, Li (bcc metals) and solid Br2 (Cmca, the low-temperature crystal) at the hull settings
(80 Ry, PBE-D3 zero), plus one fixed-occupation SCF per candidate on its converged geometry, which prints the PBE gap.
Run from repo root: discovery/.venv/bin/python -I overnight/make_window.py"""
import re, sys
from pathlib import Path
from pymatgen.core import Lattice, Structure
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inputs as mi

mi.OUT = Path(__file__).resolve().parent / "m5_window/jobs"
ON = mi.OUT.parents[1]
# metals: denser k (0.10 1/A) since the fixed 0.005 Ry smearing is small; Br2: experimental cell, Wyckoff 8f (Powell 1984)
mi.write("01_Na_bcc_ref_vcrelax", Structure.from_spacegroup("Im-3m", Lattice.cubic(4.23), ["Na"], [[0, 0, 0]]), "vc-relax",
         k=mi.kmesh(Structure.from_spacegroup("Im-3m", Lattice.cubic(4.23), ["Na"], [[0, 0, 0]]).get_primitive_structure(), 0.10), ecut=80.0)
li = Structure.from_spacegroup("Im-3m", Lattice.cubic(3.44), ["Li"], [[0, 0, 0]])
mi.write("02_Li_bcc_ref_vcrelax", li, "vc-relax", k=mi.kmesh(li.get_primitive_structure(), 0.10), ecut=80.0)
br2 = Structure.from_spacegroup(64, Lattice.orthorhombic(6.67, 4.48, 8.72), ["Br"], [[0, 0.1015, 0.1328]])
assert len(br2) == 8, len(br2)
mi.write("03_Br2_Cmca_ref_vcrelax", br2, "vc-relax", ecut=80.0)
for i, (lab, out) in enumerate((("NaSbBr4", "night3/jobs/09_NaSbBr4_TP01_hit_vcrelax.out"),
                                 ("LiSbBr4", "night3/jobs/02_LiSbBr4_TP01_hit_vcrelax.out")), 4):
    txt = (ON / out).read_text(errors="replace")
    assert "End final coordinates" in txt, out
    nelec = float(re.search(r"number of electrons\s*=\s*([\d.]+)", txt).group(1))
    name = f"0{i}_{lab}_TP01_gap_scf"
    mi.write(name, mi.final_structure(ON / out), "scf", ecut=80.0)
    p = mi.OUT / f"{name}.in"
    t = p.read_text()
    old = "  occupations='smearing'\n  smearing='gaussian'\n  degauss=0.005\n"
    assert t.count(old) == 1
    p.write_text(t.replace(old, f"  occupations='fixed'\n  nbnd={int(nelec) // 2 + 16}\n"))
