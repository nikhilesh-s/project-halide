"""M5 batch m5_oct17c (2026-10-17): the gaps the per-material reviews found (discovery/results/review_oct17/), best value first.
Run once from repo root: discovery/.venv/bin/python -I overnight/make_oct17c.py   (idempotent)
1. Moves slow Air jobs that block a gate to the M5 (the Air drops to Light): Na3BiCl6 rerun, second listed NaSbCl4, NaSbBr4
   100 Ry cutoff SCFs, Na5BiCl8 rerun, GaI3 (NaGaI4 G3 blocker, with its clean partial segments) and NaGaI3.
2. Uncomputed listed polymorphs / references: NaBiBr4 agm003733888, NaGaI4 agm077817710, YBr3 R-3 agm002139115.
3. D3-BJ / no-D3 re-relaxations for NaBiBr4 (TP + Na3BiBr6; BiBr3 is in m5_oct20, NaBr in m5_oct17b).
4. 100 Ry cutoff SCFs (G4 only varied k) for LiSbBr4, LiBiBr4, NaBiBr4 and their lowest references."""
import re, subprocess, sys
from pathlib import Path
O = Path(__file__).resolve().parent; ROOT = O.parent
sys.path.insert(0, str(O))
import make_inputs as mi
from make_nights6plus import start
from reprioritize_oct16 import relax, cutoff_scf, have, MODES

BATCH = "m5_oct17c"
MOVE = ["air/jobs/22_Na3BiCl6_agm016680864_vcrelax", "air/jobs/24_NaSbCl4_agm003736684_listed2_vcrelax",
        "air/jobs/25_NaSbBr4_TP01_hit_100Ry_scf", "air/jobs/26_NaSbBr4_agm003733885_best_existing_100Ry_scf",
        "air/jobs/27_NaBr_agm003157587_ref_100Ry_scf", "air/jobs/28_SbBr3_agm002139075_ref_100Ry_scf",
        "air_heavy/jobs/05_Na5BiCl8_agm006092210_vcrelax", "later/jobs/09_GaI3_agm003220383_ref_vcrelax",
        "later/jobs/10_NaGaI3_agm005838860_ternary_vcrelax",
        "later/jobs/07_NaGaI4_agm005755304_alex_hull0_vcrelax"]  # stale-restart output (review hardening 2026-10-17): rerun clean
NEW = []  # parked 2026-10-17 (top-5 plan; were NaBiBr4 listed3, NaGaI4 listed, YBr3 R-3): overnight/parked/m5_oct17c
_PARKED_NEW = [("NaBiBr4_agm003733888_listed3", 80), ("NaGaI4_agm077817710_listed", 80), ("YBr3_Rm3_agm002139115", 80)]
D3SET = []  # parked 2026-10-17 (top-5 plan)
_PARKED_D3SET = [("NaBiBr4_TP01_CsAlCl4type_hit", "m5_oct16/jobs/00a_NaBiBr4_TP01_CsAlCl4type_hit_vcrelax"),
         ("Na3BiBr6_agm053446239", "m5_oct16/jobs/00d_Na3BiBr6_agm053446239_ternary_vcrelax")]
# 2026-10-17 top-5 plan: D3-BJ for the remaining Li-bromide competitors (m5_oct20 covered only TP + binaries)
BJ_ONLY = [("Li3SbBr6_agm006027634", "air_heavy/jobs/02_Li3SbBr6_agm006027634_ternary_vcrelax"),
           ("LiSbBr4_agm005739369_best_existing", "air/jobs/10_LiSbBr4_agm005739369_best_existing_vcrelax"),
           ("Li3BiBr6_agm005961056", "air_heavy/jobs/03_Li3BiBr6_agm005961056_ternary_vcrelax"),
           ("LiBiBr4_agm005739442_best_existing", "air/jobs/13_LiBiBr4_agm005739442_best_existing_vcrelax")]
CUTOFF = ["night3/jobs/02_LiSbBr4_TP01_hit_vcrelax", "night3/jobs/06_LiBiBr4_TP01_hit_vcrelax",
          "m5_oct15/jobs/10l_LiBr_agm002091979_ref_vcrelax",
          "air/jobs/14_BiBr3_agm002138960_ref_vcrelax"]  # SbBr3 and NaBr 100 Ry come from the moved air/27-28

if __name__ == "__main__":
    d = O / BATCH / "jobs"; d.mkdir(parents=True, exist_ok=True)
    n = max([int(m.group(1)) for f in d.glob("*.in") if (m := re.match(r"(\d+)_", f.name))] or [0])
    for src in MOVE:
        inp = O / (src + ".in")
        if not inp.exists(): continue
        n += 1; stem = f"{n:02d}_" + re.sub(r"^\d+_", "", inp.stem)
        for f in sorted(inp.parent.glob(inp.stem + ".*")):  # .in, .out, .out.part*
            tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(f)], cwd=ROOT, capture_output=True).returncode == 0
            dst = d / f.name.replace(inp.stem, stem, 1)
            subprocess.run(["git", "mv", str(f), str(dst)] if tracked else ["mv", str(f), str(dst)], cwd=ROOT, check=True)
        print("moved", src, "->", stem)
    for name, ecut in NEW:
        if not have(d, name + "_vcrelax"): relax(d, name, start(re.search(r"agm\d+", name).group(0)), ecut)
    for short, src in D3SET:
        s = mi.final_structure(O / (src + ".out"))
        for mode, vdw in MODES.items():
            if not have(d, f"{short}_{mode}_vcrelax"): relax(d, f"{short}_{mode}", s, 80, vdw)
    for src in CUTOFF: cutoff_scf(d, src)
    q = O / "QUEUE"
    if BATCH not in q.read_text().split(): q.write_text(q.read_text().rstrip("\n") + f"\n{BATCH}\n")
    for short, src in BJ_ONLY:
        if not have(d, f"{short}_D3BJ_vcrelax"): relax(d, f"{short}_D3BJ", mi.final_structure(O / (src + ".out")), 80, MODES["D3BJ"])
