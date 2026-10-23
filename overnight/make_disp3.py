"""Third dispersion model for the top two (2026-10-20, REFEREE_REVIEW.md M1): nonlocal vdW-DF-optB86b (input_dft='vdw-df-ob86').
Why this one: TS and MBD stop on PAW in this QE 7.5 build (PW/src/setup.f90: 'Tkatchenko-Scheffler not implemented with USPP/PAW';
MBD sets ts_vdw too). XDM would run but is still a pairwise C6/C8 add-on; a nonlocal functional changes the density functional itself
and has a chance of fixing the SbBr3 ordering (experimental Pnma above hypothetical R-3 with D3 zero, D3-BJ and no D3).
Each phase restarts from its converged PBE-D3 (zero damping) geometry; same 80/640 Ry, k-spacing 0.22, PSlibrary PAW sets;
only vdw_corr='dft-d3' is swapped for input_dft. Order: SbBr3 pair, NaSbBr4 set, LiSbBr4 set, ternaries for the full hulls,
then one D3-family gap (no-D3 listed LiSbBr4). verdict.py ignores all of these (it needs vdw_corr='dft-d3' without dftd3_version).
Run from repo root: discovery/.venv/bin/python -I overnight/make_disp3.py   (idempotent)"""
import sys
from pathlib import Path
O = Path(__file__).resolve().parent
sys.path.insert(0, str(O))
import make_inputs as mi

BATCH = "m5_disp3"
NL = "  input_dft='vdw-df-ob86'\n"
JOBS = [  # (short name, base PBE-D3 vc-relax, mode)
    ("SbBr3_agm002139075", "air/jobs/12_SbBr3_agm002139075_ref_vcrelax", "OB86"),  # R-3, lowest SbBr3 in every hull so far
    ("SbBr3_Pnma_expt_agm003218591", "m5_oct15/jobs/10_SbBr3_Pnma_expt_agm003218591_vcrelax", "OB86"),  # experimental form
    ("NaSbBr4_TP01_hit", "night3/jobs/09_NaSbBr4_TP01_hit_vcrelax", "OB86"),
    ("NaBr_agm003157587", "air/jobs/15_NaBr_agm003157587_ref_vcrelax", "OB86"),
    ("NaSbBr4_agm003733885_best_existing", "night3/jobs/10_NaSbBr4_agm003733885_best_existing_vcrelax", "OB86"),
    ("LiSbBr4_TP01_hit", "night3/jobs/02_LiSbBr4_TP01_hit_vcrelax", "OB86"),
    ("LiBr_wurtzite_agm002091979", "m5_oct15/jobs/10l_LiBr_agm002091979_ref_vcrelax", "OB86"),
    ("LiBr_rocksalt_agm003157564", "m5_oct15/jobs/02_LiBr_rocksalt_agm003157564_vcrelax", "OB86"),
    ("LiSbBr4_agm005739369_best_existing", "air/jobs/10_LiSbBr4_agm005739369_best_existing_vcrelax", "OB86"),
    ("Na3SbBr6_agm016963778", "night3/jobs/12_Na3SbBr6_agm016963778_vcrelax", "OB86"),
    ("Na2SbBr5_agm005863603", "m5_oct15/jobs/10j_Na2SbBr5_agm005863603_ternary_vcrelax", "OB86"),
    ("Li3SbBr6_agm006027634", "air_heavy/jobs/02_Li3SbBr6_agm006027634_ternary_vcrelax", "OB86"),
    ("LiSbBr4_agm005739369_best_existing", "air/jobs/10_LiSbBr4_agm005739369_best_existing_vcrelax", "noD3"),  # D3BJ in m5_top2/16
]

if __name__ == "__main__":
    mi.OUT = O / BATCH / "jobs"
    mi.OUT.mkdir(parents=True, exist_ok=True)
    for n, (short, src, mode) in enumerate(JOBS, 1):
        name = f"{n:02d}_{short}_{mode}_vcrelax"
        if (mi.OUT / f"{name}.in").exists():
            continue
        out = O / (src + ".out")
        assert "JOB DONE" in out.read_text(errors="replace") and "bfgs converged" in out.read_text(errors="replace"), src
        mi.write(name, mi.final_structure(out), "vc-relax", ecut=80.0)
        p = mi.OUT / f"{name}.in"
        t = p.read_text()
        assert t.count("  vdw_corr='dft-d3'\n") == 1
        p.write_text(t.replace("  vdw_corr='dft-d3'\n", NL if mode == "OB86" else ""))
    q = O / "QUEUE"
    if BATCH not in q.read_text().split():
        q.write_text(f"{BATCH}\n" + q.read_text())
