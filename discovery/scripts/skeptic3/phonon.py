"""Attack 4: phonopy finite displacement (0.01 A), 2x2x2 supercell (96 atoms; force evaluations only), MACE-MPA-0 + D3 float64,
on the rattled+tightly relaxed hits. Frequencies at the 8 commensurate q-points only. If min < -0.3 THz: follow that mode once
(real displacement in the smallest commensurate cell, <= 24 atoms, max atomic displacement 0.15 A, both signs), relax, report.
Run from repo root after run.py."""
import json, sys, time, itertools, numpy as np
sys.path.insert(0, "discovery/scripts/skeptic3")
from pathlib import Path
from multiprocessing import Pool
from run import calc, job, OUT, HITS
PH = Path("discovery/results/skeptic3/phonon"); PH.mkdir(parents=True, exist_ok=True)
QS = [np.array(q)/2 for q in itertools.product((0, 1), repeat=3)]

def phonon(h):
    f = PH/f"{h}.json"
    if f.exists(): return json.load(open(f))
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from pymatgen.core import Structure
    from ase import Atoms
    s = Structure.from_dict(json.load(open(OUT/f"{h}_d3.json"))["structure"]); t = time.time()
    ph = Phonopy(PhonopyAtoms(symbols=[str(x.specie) for x in s], cell=s.lattice.matrix, scaled_positions=s.frac_coords),
                 supercell_matrix=np.diag([2, 2, 2]), primitive_matrix=None)
    ph.generate_displacements(distance=0.01); F = []
    for sc in ph.supercells_with_displacements:
        a = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True); a.calc = calc(True); F.append(a.get_forces())
    ph.forces = F; ph.produce_force_constants(); ph.symmetrize_force_constants()
    ph.run_qpoints(QS, with_eigenvectors=True); d = ph.get_qpoints_dict(); fr = d["frequencies"]
    iq, im = np.unravel_index(fr.argmin(), fr.shape); r = dict(id=h, n_disp=len(F), wall=time.time()-t, min_THz=float(fr.min()),
        q_min=QS[iq].tolist(), n_below_m0p3=int((fr < -0.3).sum()), per_q_min={str(q.tolist()): float(x.min()) for q, x in zip(QS, fr)},
        three_lowest_gamma=[float(x) for x in sorted(fr[0])[:4]])
    if fr.min() < -0.3:                                   # follow once
        q = QS[iq]; ev = d["eigenvectors"][iq][:, im].reshape(-1, 3); S = [i for i in range(3) if q[i] > 0]
        M = np.eye(3, dtype=int)                          # smallest cell in which q is commensurate (det 1 or 2)
        if S:
            M[S[0], S[0]] = 2
            for j in S[1:]: M[j, S[0]] = 1                # lattice vector e_j + e_i
        big = s.copy(); big.make_supercell(M); mass = np.array([x.specie.atomic_mass for x in s]); u = []
        for x in big:                                     # match each supercell atom to its primitive atom (by fractional position)
            fp = s.lattice.get_fractional_coords(x.coords); k = int(np.argmin([np.abs((fp-y.frac_coords+0.5) % 1-0.5).sum() + 9*(y.specie != x.specie) for y in s]))
            u.append(ev[k]*np.exp(2j*np.pi*q @ fp)/np.sqrt(mass[k]))
        u = np.array(u); ph0 = max(np.linspace(0, np.pi, 180, endpoint=False), key=lambda p: np.linalg.norm((u*np.exp(1j*p)).real))
        u = (u*np.exp(1j*ph0)).real; u *= 0.15/np.abs(u).max(); r["follow"] = {}
        E0 = json.load(open(OUT/f"{h}_d3.json"))["E"]
        for sign in (1, -1):
            b = big.copy()
            for i in range(len(b)): b.translate_sites([i], sign*u[i], frac_coords=False)
            j = job((f"{h}_follow{'p' if sign > 0 else 'm'}_d3", b.as_dict(), True, True, 600))
            r["follow"]["p" if sign > 0 else "m"] = dict(n=j["n"], dE_meV=1000*(j["E"]-E0), conv=j["conv"], V=j["V"])
    json.dump(r, open(f, "w")); print(r, flush=True); return r

if __name__ == "__main__":
    with Pool(4) as p: p.map(phonon, sys.argv[1:] or [h for v in HITS.values() for h in v])   # argv: also run baselines, e.g. agm084661913
