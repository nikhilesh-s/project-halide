"""S0 gate. Usage: python check_s0.py $RUN_DIR/S0/scf_TP-NaBiCl4-01.out   (exit 0 = PASS)"""
import json, re, sys
from pathlib import Path
exp = json.loads((Path(__file__).parent / "expected.json").read_text())
txt = Path(sys.argv[1]).read_text()
assert "JOB DONE" in txt, "pw.x did not finish"
e = float(re.findall(r"^!\s+total energy\s+=\s+(\S+) Ry", txt, re.M)[-1])
d = (e - exp["laptop_final_energy_Ry"]) * 13605.693123 / exp["natoms"]
ok = abs(d) < exp["tolerance_meV_per_atom"]
print(f"cluster {e:.8f} Ry  laptop {exp['laptop_final_energy_Ry']:.8f} Ry  diff {d:+.4f} meV/atom  ->", "PASS" if ok else "FAIL: STOP, fix the environment")
sys.exit(0 if ok else 1)
