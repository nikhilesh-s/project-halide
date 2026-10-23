"""Regenerate every input in anish_packet/ from the relaxed QE outputs in this repo.
Run from repo root (laptop, seconds, single thread):  discovery/.venv/bin/python -I anish_packet/make_packet.py
Reads only; writes only under anish_packet/. QE input format is reused from overnight/make_inputs.py (same
pseudopotentials, 60/480 Ry, PBE-D3, Gaussian smearing 0.005 Ry)."""
import gzip, hashlib, json, re, shutil, types
from pathlib import Path
import numpy as np
import phonopy
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROOT = Path(__file__).resolve().parents[1]
PK = ROOT / "anish_packet"
# overnight/make_inputs.py with Br added to its hard-coded species tuple, exactly as discovery/scripts/gate3/make_night3.py
# wrote the night-3 bromide inputs (Cl-only inputs are unchanged by this).
mi = types.ModuleType("make_inputs"); mi.__file__ = str(ROOT / "overnight/make_inputs.py")
exec(compile((ROOT / "overnight/make_inputs.py").read_text().replace('("Li", "Na", "Cl", "Bi", "Sb")', '("Li", "Na", "Cl", "Br", "Bi", "Sb")'),
             mi.__file__, "exec"), mi.__dict__)
mi.MASS["Br"] = 79.904; mi.UPF["Br"] = "Br.pbe-n-kjpaw_psl.1.0.0.UPF"

A = "anish_results/ionic_discovery/dft/dft_li2cabr4/"
# id -> (relaxed pw.x output, label)
SRC = {
    "TP-NaBiCl4-01": ("overnight/night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.out", "PROPOSED MODEL"),
    "agm088474086": ("overnight/night1/jobs/01_NaBiCl4_agm088474086_vcrelax.out", "Alexandria entry, re-relaxed PBE-D3"),
    "agm005743614": ("overnight/night2/jobs/00_NaBiCl4_agm005743614_rerelax_vcrelax.out", "Alexandria entry (old C2/c), re-relaxed PBE-D3"),
    "TP-NaSbCl4-01": ("overnight/night2/jobs/05_NaSbCl4_TP01_CsAlCl4type_vcrelax.out", "PROPOSED MODEL"),
    "Li2CaBr4_PBE": (A + "Li2CaBr4.out", "Alexandria agm003455102, Anish's PBE vc-relax"),
    "Li2CaBr4_PBED3": (A + "Li2CaBr4_d3.out", "Alexandria agm003455102, Anish's PBE-D3 vc-relax"),
    "TP-NaSbBr4-01": ("overnight/night3/jobs/09_NaSbBr4_TP01_hit_vcrelax.out", "PROPOSED MODEL"),      # added 2026-10-14 (S4), 80 Ry
    "TP-LiSbBr4-01": ("overnight/night3/jobs/02_LiSbBr4_TP01_hit_vcrelax.out", "PROPOSED MODEL"),      # added 2026-10-16 (S5), 80 Ry
    "TP-LiBiBr4-01": ("overnight/night3/jobs/06_LiBiBr4_TP01_hit_vcrelax.out", "PROPOSED MODEL"),      # added 2026-10-16 (S5), 80 Ry
    "TP-NaBiBr4-01": ("overnight/m5_oct16/jobs/00a_NaBiBr4_TP01_CsAlCl4type_hit_vcrelax.out", "PROPOSED MODEL"),  # S5 optional, 80 Ry
}
DISP = 0.02          # A, finite displacement
MIN_LEN = 11.0       # A, minimum supercell lattice length for phonons
STEPS = {300: 75000, 350: 60000, 400: 50000, 500: 35000, 600: 25000}   # 2 fs steps; follows Anish's rt/fin schedules
SEEDS = [1, 2, 3, 4]
manifest = {"phonopy_version": phonopy.__version__, "displacement_A": DISP, "structures": {}, "phonons": {}, "md": {}}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def input_kmesh(out):
    lines = Path(str(out)[:-4] + ".in").read_text().splitlines()
    return [int(x) for x in lines[lines.index("K_POINTS automatic") + 1].split()[:3]]


def pw_input(folder, name, s, k, conv="1.0d-8", ecut=60.0):
    """scf input via overnight/make_inputs.write, minus max_seconds, with the requested conv_thr."""
    mi.OUT = folder; folder.mkdir(parents=True, exist_ok=True)
    mi.write(name, s, "scf", k=k, ecut=ecut)
    f = folder / f"{name}.in"
    t = f.read_text().replace("  max_seconds=__MAXSEC__\n", "").replace("conv_thr=1.0d-8", f"conv_thr={conv}")
    assert "max_seconds" not in t and f"conv_thr={conv}" in t and "tprnfor=.true." in t
    f.write_text(t)


def phonon_case(folder, sid, s, k_prim, symprec=1e-5, ecut=60.0):
    mult = [int(np.ceil(MIN_LEN / x)) for x in s.lattice.abc]
    ph = Phonopy(PhonopyAtoms(symbols=[x.specie.symbol for x in s], cell=s.lattice.matrix, scaled_positions=s.frac_coords),
                 supercell_matrix=np.diag(mult), primitive_matrix="P", symprec=symprec)
    ph.generate_displacements(distance=DISP)
    k = [int(np.ceil(a / m)) for a, m in zip(k_prim, mult)]
    cells = [ph.supercell] + list(ph.supercells_with_displacements)      # disp-000 = undisplaced reference
    for i, c in enumerate(cells):
        pw_input(folder, f"disp-{i:03d}", Structure(c.cell, c.symbols, c.scaled_positions), k, conv="1.0d-10", ecut=ecut)
    ph.save(folder / "phonopy_disp.yaml")
    sc = Structure(ph.supercell.cell, ph.supercell.symbols, ph.supercell.scaled_positions)
    info = dict(structure_id=sid, supercell=mult, natoms_supercell=len(sc), n_displacements=len(cells) - 1,
                n_scf_including_reference=len(cells), kmesh_primitive=k_prim, kmesh_supercell=k,
                supercell_lengths_A=[round(x, 3) for x in sc.lattice.abc],
                spacegroup_used=SpacegroupAnalyzer(s, 1e-3).get_space_group_symbol())
    manifest["phonons"][str(folder.relative_to(PK))] = info
    return info


def md_mult(s, w):
    """Supercell multiplicity for min_width w, by Anish's supercell_for_md rule (perpendicular widths)."""
    m = s.lattice.matrix; vol = abs(np.linalg.det(m))
    wd = np.array([vol / np.linalg.norm(np.cross(m[(i + 1) % 3], m[(i + 2) % 3])) for i in range(3)])
    return np.maximum(1, np.ceil(w / wd)).astype(int).tolist()


def md_items(path, ids, temps, mobile, width=None, steps=None, seeds=SEEDS):
    items = []; path.parent.mkdir(parents=True, exist_ok=True)
    for sid in ids:
        s = raw[sid].copy(); s.properties = {"packet_id": sid}
        w = width or 10.0           # default: smallest width (0.5 A steps) giving >= 100 atoms
        while width is None and len(s) * np.prod(md_mult(s, w)) < 100:
            w += 0.5
        mult = md_mult(s, w); n = int(len(s) * np.prod(mult))
        assert n <= 320
        d = s.as_dict()
        assert Structure.from_dict(d).properties["packet_id"] == sid
        manifest["md"].setdefault(path.name, {})[sid] = dict(min_width=w, mult=mult, natoms=n, vol_per_atom_A3=round(s.volume / len(s), 4))
        for T in temps:
            for seed in seeds:
                items.append(dict(id=sid, source="anish_packet", formula=s.composition.reduced_formula, mobile=mobile, structure=d,
                                  control=sid == "agm005743614", label=SRC[sid][1], T=T, steps=steps or STEPS[T], seed=seed,
                                  min_width=w, max_atoms=320))
    items.sort(key=lambda it: -it["steps"])      # equal-length runs share a GPU (4 consecutive items per array task)
    with open(path, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as g:
        g.write(json.dumps(items).encode())
    return items


# ---- structures -----------------------------------------------------------------------------------------------
(PK / "structures").mkdir(parents=True, exist_ok=True)
raw, sym = {}, {}
for sid, (out, label) in SRC.items():
    s = raw[sid] = mi.final_structure(ROOT / out)
    CifWriter(s).write_file(PK / f"structures/{sid}_relaxed.cif")
    rec = dict(label=label, source_output=out, source_output_sha256=sha(ROOT / out), natoms=len(s), formula=s.composition.formula,
               vol_per_atom_A3=round(s.volume / len(s), 4), spacegroup_strict=SpacegroupAnalyzer(s, 1e-5).get_space_group_symbol(),
               spacegroup_loose_0p05=SpacegroupAnalyzer(s, 0.05).get_space_group_symbol(),
               cif=f"structures/{sid}_relaxed.cif", cif_sha256=sha(PK / f"structures/{sid}_relaxed.cif"))
    if sid.startswith("TP-"):      # relaxed TP cells are only approximately symmetric: phonons use the symmetrised cell
        sga = SpacegroupAnalyzer(s, 0.05)
        r = sga.get_refined_structure()
        if len(r) != len(s):        # centred group (TP-NaSbBr4-01, C2/c): keep the 12-atom primitive setting
            r = sga.get_primitive_standard_structure()
        sym[sid] = r
        assert len(r) == len(s)
        CifWriter(r, symprec=1e-3).write_file(PK / f"structures/{sid}_symmetrized.cif")
        from pymatgen.analysis.structure_matcher import StructureMatcher
        rms = StructureMatcher(primitive_cell=False, scale=False).get_rms_dist(s, r)
        rec.update(symmetrized_cif=f"structures/{sid}_symmetrized.cif", symmetrized_cif_sha256=sha(PK / f"structures/{sid}_symmetrized.cif"),
                   symmetrization_rms_and_max_shift_normalised=[round(float(x), 5) for x in rms],
                   symmetrization_max_shift_A=round(float(rms[1] * (s.volume / len(s)) ** (1 / 3)), 4),
                   symmetrized_abc=[round(x, 4) for x in r.lattice.abc], symmetrized_angles=[round(x, 3) for x in r.lattice.angles],
                   relaxed_abc=[round(x, 4) for x in s.lattice.abc], relaxed_angles=[round(x, 3) for x in s.lattice.angles])
    manifest["structures"][sid] = rec

# ---- pseudopotentials (identical files to the laptop runs) -------------------------------------------------------
(PK / "pseudo").mkdir(exist_ok=True)
md5 = dict(l.split()[::-1] for l in (ROOT / "overnight/pseudo_md5.txt").read_text().splitlines())
for el in ("Na", "Cl", "Bi", "Sb", "Br", "Li"):
    f = ROOT / "overnight/pseudo" / mi.UPF[el]
    assert hashlib.md5(f.read_bytes()).hexdigest() == md5[f.name]
    shutil.copy2(f, PK / "pseudo" / f.name)
(PK / "pseudo/MD5SUMS").write_text("".join(f"{md5[mi.UPF[e]]}  {mi.UPF[e]}\n" for e in ("Na", "Cl", "Bi", "Sb", "Br", "Li")))

# ---- S0: control SCF at the relaxed geometry ---------------------------------------------------------------------
tp_out = ROOT / SRC["TP-NaBiCl4-01"][0]
pw_input(PK / "S0", "scf_TP-NaBiCl4-01", raw["TP-NaBiCl4-01"], input_kmesh(tp_out))
e_ref = json.loads((ROOT / "overnight/night2/summary.json").read_text())["01_NaBiCl4_TP01_KBiF4type_vcrelax"]["final_energy_Ry"]
assert f"{e_ref:.8f}" in tp_out.read_text()
(PK / "S0/expected.json").write_text(json.dumps(dict(input="scf_TP-NaBiCl4-01.in", natoms=12, laptop_final_energy_Ry=e_ref,
    tolerance_meV_per_atom=0.1, source="overnight/night2/summary.json"), indent=1))

# ---- S1 / S3 phonons ---------------------------------------------------------------------------------------------
k_tp, k_agm = input_kmesh(tp_out), input_kmesh(ROOT / SRC["agm088474086"][0])
assert mi.kmesh(sym["TP-NaBiCl4-01"]) == k_tp, "symmetrised cell changes the k-mesh"
phonon_case(PK / "S1/TP-NaBiCl4-01", "TP-NaBiCl4-01 (symmetrised)", sym["TP-NaBiCl4-01"], k_tp)
phonon_case(PK / "S1/TP-NaBiCl4-01_P1fallback", "TP-NaBiCl4-01 (as relaxed, no symmetry)", raw["TP-NaBiCl4-01"], k_tp)
phonon_case(PK / "S1/agm088474086", "agm088474086", raw["agm088474086"], k_agm)
phonon_case(PK / "S3/TP-NaSbCl4-01", "TP-NaSbCl4-01 (symmetrised)", sym["TP-NaSbCl4-01"], mi.kmesh(sym["TP-NaSbCl4-01"]))

# ---- S2 / S3 MD items (format of candidates/rt_items.json.gz) ------------------------------------------------------
T5 = [300, 350, 400, 500, 600]
md_items(PK / "S2/md_items.json.gz", ["TP-NaBiCl4-01", "agm088474086", "agm005743614"], T5, "Na")
md_items(PK / "S2/smoke_items.json.gz", ["TP-NaBiCl4-01"], [300], "Na", steps=500, seeds=[1])
md_items(PK / "S3/md_nasbcl4_items.json.gz", ["TP-NaSbCl4-01"], T5, "Na")
md_items(PK / "S3/md_li2cabr4_items.json.gz", ["Li2CaBr4_PBE", "Li2CaBr4_PBED3"], [300, 350], "Li", width=12.0)

# ---- S4 (added 2026-10-14): extra stages, run after S2 and only if time remains -----------------------------------------
S4 = PK / "S4"
MAXSEC = {"relax": 20000, "rival": 27000}      # pw.x stops itself (and still writes JOB DONE) before the --time limit
manifest["s4"] = {}


def verbatim(folder, name, src_in, maxsec):
    """Copy a laptop input unchanged except max_seconds (run_pw.sh fills __PSEUDO__ and __TMP__)."""
    folder.mkdir(parents=True, exist_ok=True)
    t = (ROOT / src_in).read_text()
    assert t.count("max_seconds=__MAXSEC__") == 1
    (folder / f"{name}.in").write_text(t.replace("__MAXSEC__", str(maxsec)))
    manifest["s4"][f"{folder.name}/{name}.in"] = dict(source_input=src_in, source_input_sha256=sha(ROOT / src_in))


def tight_relax(folder, name, s, k, nosym):
    """vc-relax with forces 10x, energy 10x, pressure 5x and SCF 100x tighter than the laptop relaxations."""
    pw_input(folder, name, s, k, conv="1.0d-10")
    f = folder / f"{name}.in"
    t = f.read_text()
    for a, b in [("calculation='scf'", "calculation='vc-relax'"), ("forc_conv_thr=1.0d-3", "forc_conv_thr=1.0d-4"),
                 ("etot_conv_thr=1.0d-5", "etot_conv_thr=1.0d-6"), ("press_conv_thr=0.5", "press_conv_thr=0.1"),
                 ("nstep=150\n", f"nstep=150\n  max_seconds={MAXSEC['relax']}\n")] + ([("ibrav=0\n", "ibrav=0\n  nosym=.true.\n")] if nosym else []):
        assert t.count(a) == 1, a
        t = t.replace(a, b)
    f.write_text(t)
    manifest["s4"][f"{folder.name}/{name}.in"] = dict(structure=name.split("_", 1)[1], nosym=nosym, kmesh=k)


# C9: the laptop relaxation of TP-NaBiCl4-01 repeated on the cluster from the identical input (same start geometry)
verbatim(S4 / "relax", "repeat_TP-NaBiCl4-01", "overnight/night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.in", MAXSEC["relax"])
# C10: is P2_12_12 exact? tight relax from the symmetrised cell (symmetry kept) vs from the as-relaxed P1 cell (nosym)
tight_relax(S4 / "relax", "symtight_TP-NaBiCl4-01", sym["TP-NaBiCl4-01"], k_tp, nosym=False)
tight_relax(S4 / "relax", "P1tight_TP-NaBiCl4-01", raw["TP-NaBiCl4-01"], k_tp, nosym=True)
# backup for the laptop run of the 24-atom KGaI4-type NaBiCl4 rival (skip if its result is already in git)
verbatim(S4 / "rival", "10a1_NaBiCl4_TPscan01_KGaI4type_RIVAL_vcrelax",
         "overnight/night2/jobs/10a1_NaBiCl4_TPscan01_KGaI4type_RIVAL_vcrelax.in", MAXSEC["rival"])
# phonons of TP-NaSbBr4-01 at its own 80/640 Ry setting (only if its G4 check passes)
nsb_out = ROOT / SRC["TP-NaSbBr4-01"][0]
assert "ecutwfc=80.0" in Path(str(nsb_out)[:-4] + ".in").read_text() and mi.kmesh(sym["TP-NaSbBr4-01"]) == input_kmesh(nsb_out)
phonon_case(S4 / "TP-NaSbBr4-01", "TP-NaSbBr4-01 (symmetrised)", sym["TP-NaSbBr4-01"], input_kmesh(nsb_out), ecut=80.0)
manifest["s4"]["TP-NaSbBr4-01_relaxed_energy_Ry"] = float(re.findall(r"^!\s+total energy\s+=\s+(\S+) Ry", nsb_out.read_text(), re.M)[-1])
# MACE MD (ML, not DFT): TP-NaSbBr4-01 with 8 seeds, and seeds 5-8 for the four S2/S3 Na structures; 450 and 550 K added
STEPS.update({450: 45000, 550: 30000})
T7 = sorted(STEPS)
md_items(S4 / "md_nasbbr4_items.json.gz", ["TP-NaSbBr4-01"], T7, "Na", seeds=list(range(1, 9)))
md_items(S4 / "md_more_items.json.gz", ["TP-NaBiCl4-01", "agm088474086", "agm005743614", "TP-NaSbCl4-01"], T7, "Na", seeds=[5, 6, 7, 8])

# ---- S5 (added 2026-10-16): lithium bromides, built exactly like S4/TP-NaSbBr4-01; NaBiBr4 phonons optional (go/no-go) ----
S5 = PK / "S5"
manifest["s5"] = {}
for sid in ("TP-LiSbBr4-01", "TP-LiBiBr4-01", "TP-NaBiBr4-01"):
    out = ROOT / SRC[sid][0]
    assert "ecutwfc=80.0" in Path(str(out)[:-4] + ".in").read_text() and mi.kmesh(sym[sid]) == input_kmesh(out)
    phonon_case(S5 / sid, f"{sid} (symmetrised)", sym[sid], input_kmesh(out), ecut=80.0)
    manifest["s5"][f"{sid}_relaxed_energy_Ry"] = float(re.findall(r"^!\s+total energy\s+=\s+(\S+) Ry", out.read_text(), re.M)[-1])
# MACE MD (ML, not DFT) for Li transport: S2/S3 protocol (300-600 K, seeds 1-4); only the mobile species differs
# one file per structure: LiBiBr4 gets 216 atoms by the >= 100-atom rule, so it runs 2 per GPU (NPROC=2) to stay inside 4 h
md_items(S5 / "md_lisbbr4_items.json.gz", ["TP-LiSbBr4-01"], T5, "Li")
md_items(S5 / "md_libibr4_items.json.gz", ["TP-LiBiBr4-01"], T5, "Li")

(PK / "manifest.json").write_text(json.dumps(manifest, indent=1))
print(json.dumps({k: manifest[k] for k in ("phonons", "md")}, indent=1))
for sid, r in manifest["structures"].items():
    print(sid, {k: v for k, v in r.items() if "sha" not in k and k not in ("cif", "symmetrized_cif", "source_output", "label")})
