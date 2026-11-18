"""Dispersion check for TP-LiSbBr4-01 and TP-LiBiBr4-01 (skeptic review 2026-10-16: their margin over LiBr + MBr3 exists only with D3).
Writes vc-relax inputs into overnight/m5_oct20/jobs: each phase restarted from its final PBE-D3 (zero damping) geometry and relaxed again with
(a) D3-BJ (dftd3_version=4) and (b) no dispersion (plain PBE). Same 80/640 Ry, k-spacing 0.22, pseudopotentials as the base runs.
LiBiBr4 system first (thinnest margin). Zincblende LiBr and P2_12_12_1 SbBr3 (relaxed into Pnma) are left out.
Run from repo root: discovery/.venv/bin/python -I overnight/make_dispersion_check.py"""
import sys
from pathlib import Path
O = Path(__file__).resolve().parent
sys.path.insert(0, str(O))
import make_inputs as mi

BATCH = "m5_oct20"
SOURCES = [  # (short name, base PBE-D3 vc-relax)
    ("LiBiBr4_TP01_hit", "night3/jobs/06_LiBiBr4_TP01_hit_vcrelax"),
    ("BiBr3_agm002138960", "air/jobs/14_BiBr3_agm002138960_ref_vcrelax"),
    ("BiBr3_C2m_expt_agm002330398", "m5_oct15/jobs/10i_BiBr3_C2m_expt_agm002330398_vcrelax"),
    ("LiBr_wurtzite_agm002091979", "m5_oct15/jobs/10l_LiBr_agm002091979_ref_vcrelax"),
    ("LiBr_rocksalt_agm003157564", "m5_oct15/jobs/02_LiBr_rocksalt_agm003157564_vcrelax"),
    ("LiSbBr4_TP01_hit", "night3/jobs/02_LiSbBr4_TP01_hit_vcrelax"),
    ("SbBr3_agm002139075", "air/jobs/12_SbBr3_agm002139075_ref_vcrelax"),
    ("SbBr3_agm002257875", "m5_oct15/jobs/09_SbBr3_agm002257875_vcrelax"),
    ("SbBr3_Pnma_expt_agm003218591", "m5_oct15/jobs/10_SbBr3_Pnma_expt_agm003218591_vcrelax"),
]
MODES = {"D3BJ": "  vdw_corr='dft-d3'\n  dftd3_version=4\n", "noD3": ""}

if __name__ == "__main__":
    mi.OUT = O / BATCH / "jobs"
    mi.OUT.mkdir(parents=True, exist_ok=True)
    n = 0
    for short, src in SOURCES:
        out = O / (src + ".out")
        assert "JOB DONE" in out.read_text(errors="replace"), src
        s = mi.final_structure(out)
        for mode, vdw in MODES.items():
            n += 1
            name = f"{n:02d}_{short}_{mode}_vcrelax"
            if (mi.OUT / f"{name}.in").exists():
                continue
            mi.write(name, s, "vc-relax", ecut=80.0)
            p = mi.OUT / f"{name}.in"
            t = p.read_text()
            assert t.count("  vdw_corr='dft-d3'\n") == 1
            p.write_text(t.replace("  vdw_corr='dft-d3'\n", vdw))
    q = O / "QUEUE"
    if BATCH not in q.read_text().split():
        q.write_text(q.read_text().rstrip("\n") + f"\n{BATCH}\n")
