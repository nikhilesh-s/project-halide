"""G4 convergence check for the proposed NaSbBr4 model: SCF at the relaxed 80 Ry geometries with k-spacing 0.17
(base setting 80 Ry k0.22). Inputs are the original vc-relax inputs with calculation='scf', the final geometry and a
denser k mesh, so species and every other setting stay identical. Queued first in the Air light queue (14b*).
Run: discovery/.venv/bin/python -I overnight/make_converge_nasbbr4.py"""
import re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_inputs as mi
import run_night as rn
O = mi.ROOT / "overnight"
for name, src in [("14b1_NaSbBr4_TP01", "night3/jobs/09_NaSbBr4_TP01_hit_vcrelax"),
                  ("14b2_NaSbBr4_agm003733885", "night3/jobs/10_NaSbBr4_agm003733885_best_existing_vcrelax"),
                  ("14b3_SbBr3_agm002139075", "air/jobs/12_SbBr3_agm002139075_ref_vcrelax"),
                  ("14b4_NaBr_agm003157587", "air/jobs/15_NaBr_agm003157587_ref_vcrelax"),
                  ("14b5_Na3SbBr6_agm016963778", "night3/jobs/12_Na3SbBr6_agm016963778_vcrelax")]:
    out = (O / (src + ".out")).read_text()
    assert "JOB DONE" in out and "bfgs converged" in out, src
    job = name + "_80Ry_k017_scf"
    txt = rn.with_geometry((O / (src + ".in")).read_text(), rn.last_geometry(out))
    txt = txt.replace("calculation='vc-relax'", "calculation='scf'")
    txt = re.sub(r"prefix='[^']*'", f"prefix='{job}'", txt)
    txt = re.sub(r"outdir='__TMP__/[^']*'", f"outdir='__TMP__/{job}'", txt)
    k = mi.kmesh(mi.final_structure(O / (src + ".out")), 0.17)
    txt = re.sub(r"(K_POINTS automatic\n)\s*\d+ \d+ \d+", lambda m: m.group(1) + "  " + " ".join(map(str, k)), txt)
    (O / "air/jobs" / (job + ".in")).write_text(txt)
    print(job, k)
