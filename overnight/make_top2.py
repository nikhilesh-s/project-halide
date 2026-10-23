"""2026-10-17 (the project owner): all compute goes to fully defending the top two, TP-NaSbBr4-01 and TP-LiSbBr4-01.
Run once from repo root: discovery/.venv/bin/python -I overnight/make_top2.py   (idempotent)
Gaps closed (claims_review.md, energetics_review.md, live Alexandria query 2026-10-17: every Na-Sb-Br / Li-Sb-Br entry within
~50 meV/atom of the Alexandria hull that had no DFT here):
  Air  air/        NaSbBr4 listed agm005739607 (G2), Sb metal (needed once Sb5+ phases enter), NaSbBr6, LiSbBr6
  Air  air_heavy/  Na3SbBr6 P2_1/c agm003568537 (the cell that came out lowest in the chloride analogue), LiSbBr4 listed
                   agm003740876 (G2; same SG/atom count as the hit), Na2SbBr6, Li2SbBr6
  M5   m5_top2     (first in QUEUE) NaSbBr4 D3-BJ / no-D3 set moved from m5_oct17b plus the best listed NaSbBr4 in both modes,
                   NaSbBr4 and LiSbBr4 100 Ry SCFs and Li D3-BJ competitors moved from m5_oct17c, NaSbBr5 (28 atoms) last
Displaced: air/23 Na3BiCl6 and air_heavy/04 LiBiBr4 listed2 move (with their partial outputs) to later/, which runs after."""
import re, subprocess, sys
from pathlib import Path
O = Path(__file__).resolve().parent; ROOT = O.parent
sys.path.insert(0, str(O))
import make_inputs as mi
from make_nights6plus import start
from reprioritize_oct16 import relax, have, MODES

AIR = [("NaSbBr4_agm005739607_listed2", 80), ("Sb_agm003157733_ref", 80),
       ("NaSbBr6_agm005745825_ternary", 80), ("LiSbBr6_agm005745758_ternary", 80)]
AIR_HEAVY = [("Na3SbBr6_agm003568537_ternary", 80), ("LiSbBr4_agm003740876_listed2", 80),
             ("Na2SbBr6_agm005585966_ternary", 80), ("Li2SbBr6_agm005786618_ternary", 80)]
TO_LATER = ["air/jobs/23_Na3BiCl6_agm003569174_vcrelax", "air_heavy/jobs/04_LiBiBr4_agm003740513_listed2_vcrelax"]
M5_MOVE = [f"m5_oct17b/jobs/{n}" for n in (
               "13_NaSbBr4_TP01_hit_D3BJ_vcrelax", "14_NaSbBr4_TP01_hit_noD3_vcrelax",
               "15_NaBr_agm003157587_D3BJ_vcrelax", "16_NaBr_agm003157587_noD3_vcrelax",
               "17_Na3SbBr6_agm016963778_D3BJ_vcrelax", "18_Na3SbBr6_agm016963778_noD3_vcrelax",
               "19_Na2SbBr5_agm005863603_D3BJ_vcrelax", "20_Na2SbBr5_agm005863603_noD3_vcrelax")] + \
          [f"m5_oct17c/jobs/{n}" for n in (
               "03_NaSbBr4_TP01_hit_100Ry_scf", "04_NaSbBr4_agm003733885_best_existing_100Ry_scf",
               "05_NaBr_agm003157587_ref_100Ry_scf", "06_SbBr3_agm002139075_ref_100Ry_scf",
               "17_LiSbBr4_TP01_hit_100Ry_scf", "20_LiBr_agm002091979_ref_100Ry_scf",
               "22_Li3SbBr6_agm006027634_D3BJ_vcrelax", "23_LiSbBr4_agm005739369_best_existing_D3BJ_vcrelax")]
BEST_NASBBR4 = "night3/jobs/10_NaSbBr4_agm003733885_best_existing_vcrelax"
M5_LAST = [("NaSbBr5_agm005795742_ternary", 80)]


def move(src, d):
    """git mv a job (.in, .out, .out.part*) into d with the next free number; resumable outputs travel with it."""
    inp = O / (src + ".in")
    if not inp.exists(): return
    stem = f"{mi_next(d):02d}_" + re.sub(r"^\d+_", "", inp.stem)
    for f in sorted(inp.parent.glob(inp.stem + ".*")):
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(f)], cwd=ROOT, capture_output=True).returncode == 0
        dst = d / f.name.replace(inp.stem, stem, 1)
        subprocess.run(["git", "mv", str(f), str(dst)] if tracked else ["mv", str(f), str(dst)], cwd=ROOT, check=True)
    print("moved", src, "->", d.parent.name, stem)


def mi_next(d):
    return max([int(m.group(1)) for f in d.glob("*.in") if (m := re.match(r"(\d+)_", f.name))] or [0]) + 1


def fresh(d, jobs):
    for name, ecut in jobs:
        if not have(d, name + "_vcrelax"): relax(d, name, start(re.search(r"agm\d+", name).group(0)), ecut); print("new", d.parent.name, name)


if __name__ == "__main__":
    jd = lambda q: (O / q / "jobs")
    for src in TO_LATER: move(src, jd("later"))
    fresh(jd("air"), AIR)
    fresh(jd("air_heavy"), AIR_HEAVY)
    m5 = jd("m5_top2"); m5.mkdir(parents=True, exist_ok=True)
    for src in M5_MOVE: move(src, m5)
    s = mi.final_structure(O / (BEST_NASBBR4 + ".out"))
    for mode, vdw in MODES.items():
        if not have(m5, f"NaSbBr4_agm003733885_best_existing_{mode}_vcrelax"): relax(m5, f"NaSbBr4_agm003733885_best_existing_{mode}", s, 80, vdw)
    fresh(m5, M5_LAST)
    q = O / "QUEUE"; names = q.read_text().split()
    if names[0] != "m5_top2": q.write_text("\n".join(["m5_top2"] + [n for n in names if n != "m5_top2"]) + "\n")
