"""Reprioritize all queues by candidate strength (2026-10-16, after the skeptic review in discovery/results/skeptic_oct16/README.md).
Run once from repo root: discovery/.venv/bin/python -I overnight/reprioritize_oct16.py   (idempotent)

Strength order: NaSbBr4, NaSbCl4, NaBiCl4 (hold) > LiSbBr4, LiBiBr4 (D3-dependent) > NaBiBr4 (next to count) > LiYBr4, LiScBr4, NaGaI4.
1. Parks jobs that serve dead leads (LiInCl4, LiSbCl4) or competitors already shown irrelevant (40-atom Li3MBr6 sit ~40 meV/atom
   above the join): moved to overnight/parked/<batch>/, not deleted.
2. M5 QUEUE: m5_oct16 (strong-candidate G4 SCFs, wurtzite LiBr, NaBiBr4 G4) -> m5_oct20 (Li dispersion check) -> m5_oct19 (LiYBr4,
   NaGaI4, ScBr3).
3. Air, stream A (`air` queue): competitor gaps of the three strong candidates. Stream B (`air_heavy`, then `later`): LiBiBr4
   second listed polymorph, Na5BiCl8, NaBiCl4 vs KGaI4-rival dispersion check. The `later` jobs (GaI3, NaGaI3, Li3ScBr6) run after."""
import re, subprocess, sys
from pathlib import Path
O = Path(__file__).resolve().parent; ROOT = O.parent
sys.path.insert(0, str(O))
import make_inputs as mi
import run_night as rn
from make_nights6plus import start
from make_dispersion_check import MODES

PARK = {"m5_oct16": ["40_LiInCl4_", "51_Li3SbBr6_P4332", "52_InCl3_", "54_Li3BiBr6_P4332"],
        "m5_oct19": ["03_LiSbCl4_", "04_LiSbCl4_", "05_SbCl3_agm002142958_80Ry", "06_Li2SbCl5_", "07_Li3SbCl6_"]}
QUEUE = ["night2", "m5_oct15", "m5_oct16", "m5_oct17", "m5_oct18", "m5_oct20", "m5_oct19"]

# stream A: (job name, Alexandria ID, ecut) vc-relax from the Alexandria cell
AIR_A = [("Na3BiCl6_agm016680864", 60),        # 10 at: Alexandria-lowest Na3BiCl6, never computed; NaBiCl4's tightest competitor (-12.8)
         ("Na3BiCl6_agm003569174", 60),        # 20 at: second Na3BiCl6 cell, never computed
         ("NaSbCl4_agm003736684_listed2", 60)]  # 12 at: second listed NaSbCl4, only ranked by MACE (skeptic finding 6)
# stream A: 100/800 Ry SCF at the final geometry, same k-mesh: NaSbBr4 cutoff check (G4 only varied k)
CUTOFF = ["night3/jobs/09_NaSbBr4_TP01_hit_vcrelax", "night3/jobs/10_NaSbBr4_agm003733885_best_existing_vcrelax",
          "air/jobs/15_NaBr_agm003157587_ref_vcrelax", "air/jobs/12_SbBr3_agm002139075_ref_vcrelax"]
# stream B
AIR_B = [("LiBiBr4_agm003740513_listed2", 80),  # 12 at: second listed LiBiBr4, never computed (skeptic finding 5)
         ("Na5BiCl8_agm006092210", 60)]          # 14 at: Alexandria +55 meV/atom, never computed
# stream A, after the above: NaSbCl4 dispersion check (audit 2026-10-16: -6.29 -> +5.30 with D3 removed at frozen geometry)
NASBCL4_D3 = [("NaSbCl4_TP01_CsAlCl4type", "night2/jobs/05_NaSbCl4_TP01_CsAlCl4type_vcrelax"),
              ("NaCl_agm003157609", "night2/jobs/09_NaCl_agm003157609_vcrelax"),
              ("SbCl3_agm002205242", "night2/jobs/07_SbCl3_agm002205242_vcrelax"),
              ("SbCl3_agm002142958", "night2/jobs/08_SbCl3_agm002142958_vcrelax"),
              ("SbCl3_Pnma_expt_agm003219609", "m5_oct15/jobs/05_SbCl3_Pnma_expt_agm003219609_vcrelax"),
              ("Na3SbCl6_agm003569162", "m5_oct15/jobs/11_Na3SbCl6_agm003569162_vcrelax")]
# 2026-10-17 backlog: real D3-BJ / no-D3 re-relaxations for the two candidates that held under the frozen-geometry estimate.
# Shared references already queued elsewhere are not repeated (SbBr3 x3 in m5_oct20, NaCl in air/31-32, TP-NaBiCl4-01 in air_heavy/06-07).
NASBBR4_D3 = [("NaSbBr4_TP01_hit", "night3/jobs/09_NaSbBr4_TP01_hit_vcrelax"),              # stream A, 80 Ry
              ("NaBr_agm003157587", "air/jobs/15_NaBr_agm003157587_ref_vcrelax"),
              ("Na3SbBr6_agm016963778", "night3/jobs/12_Na3SbBr6_agm016963778_vcrelax"),
              ("Na2SbBr5_agm005863603", "m5_oct15/jobs/10j_Na2SbBr5_agm005863603_ternary_vcrelax")]
NABICL4_D3 = [("Na3BiCl6_agm067313194", "night2/jobs/02_Na3BiCl6_agm067313194_vcrelax"),    # stream B, 60 Ry
              ("BiCl3_agm002142654", "night2/jobs/10_BiCl3_agm002142654_vcrelax"),
              ("BiCl3_agm002142653", "m5_oct15/jobs/06_BiCl3_agm002142653_vcrelax"),
              ("BiCl3_Pnma_expt_agm003218787", "m5_oct15/jobs/07_BiCl3_Pnma_expt_agm003218787_vcrelax")]
RIVAL = [("NaBiCl4_TP01_KBiF4type", "night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax"),
         ("NaBiCl4_TPscan01_KGaI4type_RIVAL", "night2/jobs/10a1_NaBiCl4_TPscan01_KGaI4type_RIVAL_vcrelax")]


def nxt(d):
    return max([int(m.group(1)) for f in d.glob("*.in") if (m := re.match(r"(\d+)_", f.name))] or [0]) + 1


def have(d, stem):  # job already written (any number prefix)
    return any(re.sub(r"^\d+_", "", f.stem) == stem for f in d.glob("*.in"))


def park():
    for batch, prefixes in PARK.items():
        dst = O / "parked" / batch; dst.mkdir(parents=True, exist_ok=True)
        for f in sorted((O / batch / "jobs").glob("*")):
            if any(f.name.startswith(p) for p in prefixes):
                tracked = subprocess.run(["git", "ls-files", "--error-unmatch", str(f)], cwd=ROOT, capture_output=True).returncode == 0
                subprocess.run(["git", "mv", str(f), str(dst / f.name)] if tracked else ["mv", str(f), str(dst / f.name)], cwd=ROOT, check=True)
                print("parked", batch, f.name)


def relax(d, name, s, ecut, vdw=None):
    full = f"{nxt(d):02d}_{name}_vcrelax"
    mi.OUT = d; mi.write(full, s, "vc-relax", ecut=float(ecut))
    if vdw is not None:
        p = d / f"{full}.in"; t = p.read_text()
        assert t.count("  vdw_corr='dft-d3'\n") == 1
        p.write_text(t.replace("  vdw_corr='dft-d3'\n", vdw))


def cutoff_scf(d, src):
    out = O / (src + ".out"); assert "JOB DONE" in out.read_text(errors="replace"), src
    stem = re.sub(r"^\w+?_", "", Path(src).name, count=1).removesuffix("_vcrelax") + "_100Ry_scf"
    if have(d, stem): return
    job = f"{nxt(d):02d}_{stem}"
    txt = rn.with_geometry((O / (src + ".in")).read_text(), rn.last_geometry(out.read_text(errors="replace")))
    txt = txt.replace("calculation='vc-relax'", "calculation='scf'").replace("ecutwfc=80.0", "ecutwfc=100.0").replace("ecutrho=640.0", "ecutrho=800.0")
    txt = re.sub(r"prefix='[^']*'", f"prefix='{job}'", txt)
    txt = re.sub(r"outdir='__TMP__/[^']*'", f"outdir='__TMP__/{job}'", txt)
    assert "ecutwfc=100.0" in txt and "'scf'" in txt, src
    (d / f"{job}.in").write_text(txt); print(job)


if __name__ == "__main__":
    park()
    (O / "QUEUE").write_text("\n".join(QUEUE) + "\n")
    a = O / "air/jobs"; b = O / "air_heavy/jobs"
    for name, ecut in AIR_A:
        if not have(a, name + "_vcrelax"): relax(a, name, start(re.search(r"agm\d+", name).group(0)), ecut)
    for src in CUTOFF: cutoff_scf(a, src)
    for short, src in NASBCL4_D3:
        s = mi.final_structure(O / (src + ".out"))
        for mode, vdw in MODES.items():
            if not have(a, f"{short}_{mode}_vcrelax"): relax(a, f"{short}_{mode}", s, 60, vdw)
    for name, ecut in AIR_B:
        if not have(b, name + "_vcrelax"): relax(b, name, start(re.search(r"agm\d+", name).group(0)), ecut)
    for short, src in RIVAL:
        s = mi.final_structure(O / (src + ".out"))
        for mode, vdw in MODES.items():
            if not have(b, f"{short}_{mode}_vcrelax"): relax(b, f"{short}_{mode}", s, 60, vdw)
    for sets, d, ecut in ((NASBBR4_D3, a, 80), (NABICL4_D3, b, 60)):
        for short, src in sets:
            st = mi.final_structure(O / (src + ".out"))
            for mode, vdw in MODES.items():
                if not have(d, f"{short}_{mode}_vcrelax"): relax(d, f"{short}_{mode}", st, ecut, vdw)
