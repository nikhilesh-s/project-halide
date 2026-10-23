"""Commensurate-q phonon minimum (phonopy, 0.01 A displacements, 2x2x2 supercell, MACE-MPA-0 + D3 float64) for the top <= 5 new
contenders, best ranked first; force evaluations are spread over the 4 workers and cached; a structure whose forces are not
complete at the deadline is reported as 'not finished'.  argv: minutes"""
import sys, os, itertools; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from d3_final import contenders
PH = C/"ph"; QS = [np.array(q)/2 for q in itertools.product((0, 1), repeat=3)]
def setup(p):
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    ph = Phonopy(PhonopyAtoms(symbols=[str(x.specie) for x in p], cell=p.lattice.matrix, scaled_positions=p.frac_coords), supercell_matrix=np.diag([2, 2, 2]), primitive_matrix=None)
    ph.generate_displacements(distance=0.01); return ph
def force(job):
    import transplant_scan as T
    from ase import Atoms
    f, sym, cell, sp, dl = job
    if os.path.exists(f) or time.time() > dl: return
    a = Atoms(sym, cell=cell, scaled_positions=sp, pbc=True); a.calc = T._C; np.save(f, a.get_forces())
if __name__ == "__main__":
    PH.mkdir(exist_ok=True); dl = time.time() + float(sys.argv[1])*60
    top = [c for c in contenders() if not c["match"] and c["n"] <= 12][:5]; phs = {c["name"]: setup(c["p"]) for c in top}; jobs = []
    for n, ph in phs.items():
        jobs += [(str(PH/f"{n}_{i}.npy"), list(sc.symbols), np.array(sc.cell), np.array(sc.scaled_positions), dl) for i, sc in enumerate(ph.supercells_with_displacements)]
    print("phonon structures", [(n, len(ph.supercells_with_displacements)) for n, ph in phs.items()], "force jobs", len(jobs), flush=True)
    if jobs:
        with Pool(NPROC, initializer=_init, initargs=(True,)) as p: p.map(force, jobs, chunksize=1)
    for n, ph in phs.items():
        fs = [PH/f"{n}_{i}.npy" for i in range(len(ph.supercells_with_displacements))]
        if not all(f.exists() for f in fs): print(n, "not finished", sum(f.exists() for f in fs), "/", len(fs)); continue
        ph.forces = [np.load(f) for f in fs]; ph.produce_force_constants(); ph.symmetrize_force_constants(); ph.run_qpoints(QS); fr = ph.get_qpoints_dict()["frequencies"]
        r = dict(name=n, n_disp=len(fs), min_THz=float(fr.min()), q_min=QS[int(np.unravel_index(fr.argmin(), fr.shape)[0])].tolist(), n_below_m0p3=int((fr < -0.3).sum()), min_optical_gamma=float(sorted(fr[0])[3]))
        json.dump(r, open(PH/f"{n}.json", "w")); print(r)
