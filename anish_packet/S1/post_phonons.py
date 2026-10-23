"""Collect pw.x forces -> phonopy -> stability verdict.

  python post_phonons.py <packet case dir> <run case dir>
  e.g. python $PACKET/S1/post_phonons.py $PACKET/S1/agm088474086 $RUN_DIR/S1/agm088474086

Needs phonopy_disp.yaml (packet) and disp-000.out ... disp-NNN.out (run dir). disp-000 is the undisplaced supercell; its
forces are subtracted from every displaced set, so residual forces of the relaxed geometry do not leak into the force
constants. Writes into the run dir: phonon_summary.json, phonopy_params.yaml (with force constants) and, on FAIL, a
fixed-cell relax input that follows the worst unstable mode once (follow_mode_relax.in.tpl).

Rule (THz, imaginary shown negative):
  FAIL_FOLLOW_MODE  a mode below -0.3 at a q-point commensurate with the supercell other than Gamma (exact there),
                    or an optical mode below -0.3 at Gamma
  INCONCLUSIVE      below -0.3 only between commensurate points (Fourier interpolation; needs a larger supercell)
  PASS              otherwise
No non-analytical (LO-TO) correction: it shifts LO modes near Gamma only and cannot create or remove an instability.

  python post_phonons.py --selftest <packet case dir>     plumbing check with toy Lennard-Jones forces (seconds, CPU)
"""
import itertools, json, re, sys, tempfile
from pathlib import Path
import numpy as np
import phonopy

RY_AU = 13.605693123 / 0.529177210903      # Ry/bohr -> eV/A
THRESH, GAMMA_RADIUS, MESH, FOLLOW_AMP = -0.3, 0.15, [16, 16, 16], 0.15


def qe_forces(out, nat):
    txt = Path(out).read_text()
    assert "JOB DONE" in txt and "convergence NOT achieved" not in txt, f"{out}: unfinished or unconverged"
    blk = txt[txt.rindex("Forces acting on atoms"):]
    f = re.findall(r"atom\s+\d+\s+type\s+\d+\s+force\s*=\s*(\S+)\s+(\S+)\s+(\S+)", blk)[:nat]
    assert len(f) == nat, f"{out}: {len(f)} force lines, expected {nat}"
    return np.array(f, float) * RY_AU, float(re.findall(r"^!\s+total energy\s+=\s+(\S+) Ry", txt, re.M)[-1])


def analyse(case, run):
    case, run = Path(case), Path(run)
    ph = phonopy.load(case / "phonopy_disp.yaml", produce_fc=False, log_level=0)
    nat, nd = len(ph.supercell), len(ph.displacements)
    f0, e0 = qe_forces(run / "disp-000.out", nat)
    ph.forces = np.array([qe_forces(run / f"disp-{i:03d}.out", nat)[0] - f0 for i in range(1, nd + 1)])
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    mult = np.diag(ph.supercell_matrix).astype(int)

    comm = np.array(list(itertools.product(*[np.arange(m) / m for m in mult])))        # exact q-points; comm[0] is Gamma
    ph.run_qpoints(comm, with_eigenvectors=True)
    fc = ph.get_qpoints_dict()["frequencies"]
    ph.run_mesh(MESH)
    md = ph.get_mesh_dict()
    qm, fm = md["qpoints"], md["frequencies"]
    away = np.abs(qm - np.round(qm)).max(axis=1) >= GAMMA_RADIUS
    i_mesh = np.argmin(np.where(away, fm.min(axis=1), np.inf))
    ph.run_total_dos()
    dos = ph.get_total_dos_dict()
    pts = [[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [0, .5, 0], [0, 0, 0], [0, 0, .5], [.5, 0, .5], [.5, .5, .5], [0, .5, .5], [0, 0, 0]]
    ph.run_band_structure([np.linspace(a, b, 41) for a, b in zip(pts[:-1], pts[1:])])
    band = ph.get_band_structure_dict()

    gamma_optical = float(np.sort(fc[0])[3])
    i_comm = 1 + int(np.argmin(fc[1:].min(axis=1))) if len(comm) > 1 else 0
    min_comm = float(fc[i_comm].min()) if len(comm) > 1 else 0.0
    min_mesh = float(fm[i_mesh].min())
    follow = min_comm < THRESH or gamma_optical < THRESH
    verdict = "FAIL_FOLLOW_MODE" if follow else "INCONCLUSIVE" if min_mesh < THRESH else "PASS"
    summary = dict(
        case=case.name, verdict=verdict, threshold_THz=THRESH, natoms_supercell=nat, n_displacements=nd, supercell=mult.tolist(),
        reference_energy_Ry=e0, reference_energy_Ry_per_atom=e0 / nat,
        reference_max_force_eV_A=float(np.linalg.norm(f0, axis=1).max()),
        min_freq_commensurate_nonGamma_THz=min_comm, q_min_commensurate=comm[i_comm].tolist(),
        gamma_lowest6_THz=np.sort(fc[0])[:6].tolist(), gamma_lowest_optical_THz=gamma_optical,
        min_freq_mesh_away_from_Gamma_THz=min_mesh, q_min_mesh=qm[i_mesh].tolist(), mesh=MESH,
        min_freq_mesh_all_THz=float(fm.min()), max_freq_THz=float(fm.max()),
        commensurate=[dict(q=q.tolist(), lowest4_THz=np.sort(f)[:4].tolist()) for q, f in zip(comm, fc)],
        band=dict(note="generic path through the 8 zone-boundary points in reduced coordinates, not a symmetry-standard path",
                  vertices=pts, distances=[d.tolist() for d in band["distances"]], frequencies_THz=[f.tolist() for f in band["frequencies"]]),
        dos=dict(frequency_THz=dos["frequency_points"].tolist(), total=dos["total_dos"].tolist()),
        phonopy_version=phonopy.__version__)

    if follow:      # displace along the worst exact mode, max atomic shift FOLLOW_AMP A, and write a fixed-cell relax
        iq = 0 if gamma_optical < min(THRESH, min_comm) else i_comm
        ib = int(np.argsort(fc[iq])[3 if iq == 0 else 0])
        ph.run_modulations(mult.tolist(), [[comm[iq].tolist(), ib, 1.0, 0.0]])
        mods, sc = ph.get_modulations_and_supercell()
        d = np.real(mods[0]) if np.abs(np.real(mods[0])).max() > 1e-8 else np.imag(mods[0])
        d *= FOLLOW_AMP / np.linalg.norm(d, axis=1).max()
        frac = (sc.positions + d) @ np.linalg.inv(sc.cell)
        tpl = (case / "disp-000.in").read_text()
        head, tail = tpl[:tpl.index("ATOMIC_POSITIONS crystal\n")], tpl[tpl.index("K_POINTS"):]
        head = (head.replace("calculation='scf'", "calculation='relax'").replace("disp-000", "follow_mode_relax")
                .replace("conv_thr=1.0d-10", "conv_thr=1.0d-8").replace("nstep=150", "nstep=60\n  max_seconds=27000"))
        assert "calculation='relax'" in head
        pos = "".join(f"  {s} {x:.8f} {y:.8f} {z:.8f}\n" for s, (x, y, z) in zip(sc.symbols, frac))
        (run / "follow_mode_relax.in.tpl").write_text(head + "ATOMIC_POSITIONS crystal\n" + pos + tail)
        summary["follow_mode"] = dict(q=comm[iq].tolist(), band_index=ib, frequency_THz=float(fc[iq][ib]), max_shift_A=FOLLOW_AMP,
                                      input="follow_mode_relax.in.tpl", compare_with_reference_energy_Ry=e0)
    (run / "phonon_summary.json").write_text(json.dumps(summary, indent=1))
    ph.save(run / "phonopy_params.yaml", settings={"force_constants": True})
    print(json.dumps({k: v for k, v in summary.items() if k not in ("band", "dos", "commensurate")}, indent=1))
    return summary, ph


def selftest(case):
    """Toy Lennard-Jones forces written as fake pw.x outputs; force constants must equal those from the forces set directly."""
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    ph = phonopy.load(Path(case) / "phonopy_disp.yaml", produce_fc=False, log_level=0)
    cells = [ph.supercell] + list(ph.supercells_with_displacements)
    forces = []
    with tempfile.TemporaryDirectory() as tmp:
        for i, c in enumerate(cells):
            a = Atoms(c.symbols, cell=c.cell, scaled_positions=c.scaled_positions, pbc=True)
            a.calc = LennardJones(sigma=2.4, epsilon=0.02, rc=6.0)
            f = a.get_forces(); forces.append(f)
            lines = "".join(f"     atom {j+1:4d} type  1   force = {x:14.8f}{y:14.8f}{z:14.8f}\n" for j, (x, y, z) in enumerate(f / RY_AU))
            Path(tmp, f"disp-{i:03d}.out").write_text("!    total energy              =   -1.00000000 Ry\n"
                "     Forces acting on atoms (cartesian axes, Ry/au):\n\n" + lines + "\n   JOB DONE.\n")
        summary, ph1 = analyse(case, tmp)
        assert not summary.get("follow_mode") or Path(tmp, "follow_mode_relax.in.tpl").read_text().count("\n") > len(cells[0])
    ph.forces = np.array(forces[1:]) - forces[0]
    ph.produce_force_constants(); ph.symmetrize_force_constants()
    err = np.abs(ph.force_constants - ph1.force_constants).max()
    assert err < 1e-4, err
    assert len(summary["band"]["frequencies_THz"]) == 9 and len(summary["commensurate"]) == 8
    print(f"SELFTEST OK ({len(cells)} cells, {len(cells[0])} atoms, max FC difference {err:.1e} eV/A^2, toy verdict {summary['verdict']})")


if __name__ == "__main__":
    selftest(sys.argv[2]) if sys.argv[1] == "--selftest" else analyse(sys.argv[1], sys.argv[2])
