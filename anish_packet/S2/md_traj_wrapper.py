"""Run Anish's sse_md.py UNCHANGED, with two behaviours swapped in from outside:
  1. every MD trajectory is saved (the earlier runs kept only MSD numbers);
  2. the cell is held fixed at the supplied (DFT) geometry: positions are relaxed, the cell is not.
Usage is identical to sse_md.py:  python md_traj_wrapper.py --items X.json.gz --tag T --proc p --nproc n ...
Env: SSE_DIR (folder with sse_md.py, sse_common.py), TRAJ_DIR (default ./traj), TRAJ_STRIDE (default 2 -> 40 fs frames).
Each run writes TRAJ_DIR/<id>_T<temp>_s<seed>.npz: positions (frames, atoms, 3; unwrapped Cartesian A, float32),
temps (K, every sample), symbols, cell (3x3 A), frame_dt_fs, T, seed."""
import os, runpy, sys
import numpy as np

SSE_DIR = os.environ["SSE_DIR"]
sys.path.insert(0, SSE_DIR)
import sse_common as C

TRAJ_DIR = os.environ.get("TRAJ_DIR", "traj")
STRIDE = int(os.environ.get("TRAJ_STRIDE", "2"))
os.makedirs(TRAJ_DIR, exist_ok=True)
_orig_run_md = C.run_md
_cur = {"id": "unknown"}


def relax_fixed_cell(struct, calc, fmax=0.05, steps=400):
    from ase.optimize import FIRE
    from pymatgen.io.ase import AseAtomsAdaptor
    _cur["id"] = (getattr(struct, "properties", None) or {}).get("packet_id") or struct.composition.reduced_formula
    atoms = C.to_atoms(struct)
    atoms.calc = calc
    conv = FIRE(atoms, logfile=None).run(fmax=fmax, steps=steps)      # no cell filter: volume stays at the input value
    e = float(atoms.get_potential_energy())
    fm = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
    return AseAtomsAdaptor.get_structure(atoms), e, bool(conv), fm


def run_md_saving(struct, calc, T, n_steps, dt_fs=2.0, sample_every=10, seed=0):
    frames, temps, syms, vol = _orig_run_md(struct, calc, T, n_steps, dt_fs=dt_fs, sample_every=sample_every, seed=seed)
    path = os.path.join(TRAJ_DIR, f"{_cur['id']}_T{T:.0f}_s{seed}.npz")
    np.savez_compressed(path + ".tmp.npz", positions=frames[::STRIDE].astype(np.float32), temps=temps.astype(np.float32),
                        symbols=np.array(syms), cell=np.array(struct.lattice.matrix), frame_dt_fs=dt_fs * sample_every * STRIDE,
                        T=T, seed=seed)
    os.replace(path + ".tmp.npz", path)
    print(f"  saved {path} frames={len(frames[::STRIDE])} natoms={len(syms)}", flush=True)
    return frames, temps, syms, vol


C.relax = relax_fixed_cell
C.run_md = run_md_saving
runpy.run_path(os.path.join(SSE_DIR, "sse_md.py"), run_name="__main__")   # its `from sse_common import ...` now gets the two wrappers
