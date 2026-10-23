"""Write QE inputs for the overnight run from relaxed CIFs, using Anish's NaBiCl4_d3 settings verbatim.
Run from repo root: discovery/.venv/bin/python -I overnight/make_inputs.py"""
from pathlib import Path
import numpy as np
from pymatgen.core import Structure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "overnight/night1/jobs"
MASS = {"Na": 22.9898, "Cl": 35.4530, "Bi": 208.9804, "Li": 6.941, "Sb": 121.760}
MASS.update(Br=79.904, In=114.818, Sc=44.95591)  # same values as gate5/make_queue.py (night 4-5 inputs)
UPF = {"Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF", "Cl": "Cl.pbe-n-kjpaw_psl.1.0.0.UPF", "Bi": "Bi.pbe-dn-kjpaw_psl.1.0.0.UPF",
       "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF", "Sb": "Sb.pbe-n-kjpaw_psl.1.0.0.UPF"}
UPF.update(Br="Br.pbe-n-kjpaw_psl.1.0.0.UPF", In="In.pbe-dn-kjpaw_psl.1.0.0.UPF", Sc="Sc.pbe-spn-kjpaw_psl.1.0.0.UPF")
MASS["Al"], UPF["Al"] = 26.9815, "Al.pbe-n-kjpaw_psl.1.0.0.UPF"
KSPACING = 0.22  # 1/A incl. 2*pi; at least as dense on average as Anish's 4x5x4 on his NaBiCl4 cell

def kmesh(s, spacing=None):
    return [max(1, int(np.ceil(np.linalg.norm(b) / (spacing or KSPACING)))) for b in s.lattice.reciprocal_lattice.matrix]

def final_structure(out):
    """Last CELL_PARAMETERS / ATOMIC_POSITIONS block of a pw.x vc-relax output (exact cell setting, no re-reduction)."""
    import re
    txt = Path(out).read_text()
    blk = txt[txt.rindex("CELL_PARAMETERS"):]
    hdr = blk.splitlines()[0]
    lines = blk.splitlines()
    cell = np.array([[float(x) for x in l.split()] for l in lines[1:4]])
    if "alat" in hdr:
        cell *= float(re.search(r"alat=\s*([\d.]+)", hdr).group(1)) * 0.529177210903
    i = next(n for n, l in enumerate(lines) if l.startswith("ATOMIC_POSITIONS"))
    assert "crystal" in lines[i]
    sp, xyz = [], []
    for l in lines[i + 1:]:
        w = l.split()
        if len(w) != 4 or not w[0].isalpha():
            break
        sp.append(w[0]); xyz.append([float(x) for x in w[1:]])
    return Structure(cell, sp, xyz)

def write(name, s, calc, k=None, ecut=60.0):
    global OUT
    s = s.get_primitive_structure() if calc != "scf" else s
    k = k or kmesh(s)
    have = {x.symbol for x in s.composition}  # 2026-10-18: elements outside the fixed order were silently dropped (LiAlBr4 lost Al)
    els = [e for e in ("Li", "Na", "Cl", "Br", "In", "Sc", "Bi", "Sb") if e in have] + sorted(have - {"Li", "Na", "Cl", "Br", "In", "Sc", "Bi", "Sb"})
    assert all(e in MASS and e in UPF for e in els), f"no mass/pseudopotential for {[e for e in els if e not in UPF or e not in MASS]}"
    txt = f"""&CONTROL
  calculation='{calc}'
  prefix='{name}'
  pseudo_dir='__PSEUDO__'
  outdir='__TMP__/{name}'
  disk_io='low'
  tprnfor=.true.
  tstress=.true.
  forc_conv_thr=1.0d-3
  etot_conv_thr=1.0d-5
  nstep=150
  max_seconds=__MAXSEC__
/
&SYSTEM
  ibrav=0
  nat={len(s)}
  ntyp={len(els)}
  ecutwfc={ecut:.1f}
  ecutrho={8*ecut:.1f}
  vdw_corr='dft-d3'
  occupations='smearing'
  smearing='gaussian'
  degauss=0.005
/
&ELECTRONS
  conv_thr=1.0d-8
  mixing_beta=0.4
/
&IONS
/
&CELL
  press_conv_thr=0.5
/
ATOMIC_SPECIES
""" + "".join(f"  {e} {MASS[e]:.4f} {UPF[e]}\n" for e in els)
    txt += "CELL_PARAMETERS angstrom\n" + "".join("  " + " ".join(f"{x:.8f}" for x in v) + "\n" for v in s.lattice.matrix)
    txt += "ATOMIC_POSITIONS crystal\n" + "".join(f"  {site.specie.symbol} " + " ".join(f"{x:.8f}" for x in site.frac_coords) + "\n" for site in s)
    txt += f"K_POINTS automatic\n  {k[0]} {k[1]} {k[2]} 0 0 0\n"
    (OUT / f"{name}.in").write_text(txt)
    print(name, calc, s.composition.reduced_formula, "nat", len(s), "k", k, "vol/atom %.2f" % (s.volume / len(s)))

if __name__ == "__main__":
    global_out = OUT
    OUT.mkdir(parents=True, exist_ok=True)
    AR = ROOT / "anish_results/ionic_discovery/structures/qe_pbe_relaxed"
    SK = ROOT / "discovery/results/skeptic/relax"
    ref = final_structure(ROOT / "anish_results/ionic_discovery/dft/dft_nabicl4/NaBiCl4_d3.out")
    print("b norms", [round(float(np.linalg.norm(b)), 4) for b in ref.lattice.reciprocal_lattice.matrix])

    write("00_control_NaBiCl4_agm005743614_scf", ref, "scf", k=[4, 5, 4])           # must reproduce Anish's -2667.13299844 Ry
    write("01_NaBiCl4_agm088474086_vcrelax", Structure.from_file(SK / "agm088474086_d3.cif"), "vc-relax")
    write("02_NaBiCl4_agm084661913_vcrelax", Structure.from_file(SK / "agm084661913_d3.cif"), "vc-relax")
    # night 2: all PROPOSED MODELS (TP-*) plus the baselines/references needed to judge them. Same settings as night 1;
    # Li-containing cells at 80/640 Ry (Li.pbe-s PAW wants a higher cutoff; Anish used 80 Ry for Li2CaBr4).
    OUT = ROOT / "overnight/night2/jobs"; OUT.mkdir(parents=True, exist_ok=True)
    S3 = ROOT / "discovery/results/skeptic3/relax"
    q = [("01_NaBiCl4_TP01_KBiF4type", "TP-NaBiCl4-01", 60), ("02_Na3BiCl6_agm067313194", "agm067313194", 60),
         ("03_NaBiCl4_agm003736714", "agm003736714", 60), ("04_NaBiCl4_agm075995748", "agm075995748", 60),
         ("05_NaSbCl4_TP01_CsAlCl4type", "TP-NaSbCl4-01", 60), ("06_NaSbCl4_agm005743441", "agm005743441", 60),
         ("07_SbCl3_agm002205242", "agm002205242", 60), ("08_SbCl3_agm002142958", "agm002142958", 60),
         ("09_NaCl_agm003157609", "agm003157609", 60), ("10_BiCl3_agm002142654", "agm002142654", 60),
         ("11_LiBiCl4_TP01_CsAlCl4type_80Ry", "TP-LiBiCl4-01", 80), ("12_LiCl_rs_agm003157567_80Ry", "agm003157567", 80),
         ("13_LiCl_wz_agm002157414_80Ry", "agm002157414", 80), ("14_BiCl3_agm002142654_80Ry", "agm002142654", 80),
         ("15_LiBiCl4_agm003736573_80Ry", "agm003736573", 80)]
    # consistency job: re-relax the previously used structure from Anish's final geometry, so its fresh-basis
    # energy is taken at a geometry relaxed with the final basis (his run shrank 12% in volume)
    write("00_NaBiCl4_agm005743614_rerelax_vcrelax", ref, "vc-relax")
    for name, src, ecut in q:
        write(name + "_vcrelax", Structure.from_file(S3 / f"{src}_d3.cif"), "vc-relax", ecut=float(ecut))
