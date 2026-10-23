"""S6 (2026-10-18, the project owner's side): every DFT job still open for the top two, TP-NaSbBr4-01 and TP-LiSbBr4-01, moved to Anish's
cluster because the laptops are stopped. Run from repo root: discovery/.venv/bin/python -I anish_packet/make_s6.py (run once).
Writes S6/top2/disp-NNN.in (named so go.sh's existing `ph` array runs them; they are vc-relax/SCF jobs, not phonons) and
S6/top2/MANIFEST.tsv mapping each to its laptop job, so results are filed back under overnight/ by name.
Settings are the laptop inputs unchanged (QE PBE-D3, 80 Ry, PSlibrary PAW); unfinished relaxations restart from their last geometry.
Sources: forensics_oct18/NaSbBr4.md fix 1 (uncomputed SbBr3/NaBr polymorphs), make_top2.py (all other jobs)."""
import re, sys, tempfile
from pathlib import Path
R = Path(__file__).resolve().parent.parent; O = R / "overnight"
sys.path.insert(0, str(O))
import make_inputs as mi
from make_nights6plus import start
from run_night import last_geometry, with_geometry

OPEN = ["air/jobs/22_NaSbBr4_agm005739607_listed2_vcrelax", "air_heavy/jobs/04_Na3SbBr6_agm003568537_ternary_vcrelax",
        "air/jobs/23_Sb_agm003157733_ref_vcrelax", "air/jobs/24_NaSbBr6_agm005745825_ternary_vcrelax",
        "air/jobs/25_LiSbBr6_agm005745758_ternary_vcrelax", "air_heavy/jobs/05_LiSbBr4_agm003740876_listed2_vcrelax",
        "air_heavy/jobs/06_Na2SbBr6_agm005585966_ternary_vcrelax", "air_heavy/jobs/07_Li2SbBr6_agm005786618_ternary_vcrelax"] + \
       [f"m5_top2/jobs/{n}" for n in ("02_NaSbBr4_TP01_hit_noD3_vcrelax", "18_NaSbBr4_agm003733885_best_existing_noD3_vcrelax",
        "04_NaBr_agm003157587_noD3_vcrelax", "06_Na3SbBr6_agm016963778_noD3_vcrelax", "08_Na2SbBr5_agm005863603_noD3_vcrelax",
        "01_NaSbBr4_TP01_hit_D3BJ_vcrelax", "17_NaSbBr4_agm003733885_best_existing_D3BJ_vcrelax", "03_NaBr_agm003157587_D3BJ_vcrelax",
        "05_Na3SbBr6_agm016963778_D3BJ_vcrelax", "07_Na2SbBr5_agm005863603_D3BJ_vcrelax")]
NEW = ["SbBr3_agm071416924", "SbBr3_agm080347906", "SbBr3_agm071349482", "NaBr_agm080902039", "NaBr_wz_agm004462092"]
MAXSEC = 42000  # 11.7 h, inside the 12 h array time limit: QE stops cleanly and the output keeps its last geometry

if __name__ == "__main__":
    d = Path(__file__).resolve().parent / "S6" / "top2"; d.mkdir(parents=True, exist_ok=True)
    if (d / "MANIFEST.tsv").exists():  # written once: re-running after jobs finish would renumber disp-NNN under Anish
        sys.exit("S6/top2 already written; delete it by hand only if Anish has not pulled it")
    rows, n = [], 0
    for src in OPEN:
        inp, out = O / (src + ".in"), O / (src + ".out")
        t = inp.read_text()
        if out.exists():
            o = out.read_text(errors="replace")
            if "JOB DONE" in o: continue
            g = last_geometry(o)
            if g and "vc-relax" in t and "The initial density is read from file" not in o: t = with_geometry(t, g)
        rows.append((n, src, t)); n += 1
    tmp = Path(tempfile.mkdtemp())
    for name in NEW:
        mi.OUT = tmp; mi.write(f"{name}_vcrelax", start(re.search(r"agm\d+", name).group(0)), "vc-relax", ecut=80.0)
        rows.append((n, f"new/{name}_vcrelax", (tmp / f"{name}_vcrelax.in").read_text())); n += 1
    man = ["disp\tlaptop_job"]
    for i, src, t in rows:
        job = f"disp-{i:03d}"
        t = t.replace("__MAXSEC__", str(MAXSEC))
        t = re.sub(r"prefix='[^']*'", f"prefix='{job}'", t); t = re.sub(r"outdir='__TMP__/[^']*'", f"outdir='__TMP__/{job}'", t)
        assert "__PSEUDO__" in t and "__MAXSEC__" not in t, src
        (d / f"{job}.in").write_text(t); man.append(f"{job}\t{src}")
    (d / "MANIFEST.tsv").write_text("\n".join(man) + "\n")
    print(n, "jobs ->", d)
