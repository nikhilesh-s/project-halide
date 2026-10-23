"""MACE-MPA-0 screening of NaBiCl(4-x)Br(x) orderings. Resumable. Run from repo root."""
import json, time, glob, os, sys, csv
import numpy as np
from pathlib import Path
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor as AA
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.filters import FrechetCellFilter
from ase.optimize import FIRE
from mace.calculators import mace_mp

D = Path("discovery"); CAND = D/"candidates"; RES = D/"results"/"mace"; PH = D/"results"/"phonon"
RES.mkdir(parents=True, exist_ok=True); PH.mkdir(parents=True, exist_ok=True)
T0 = time.time()
calc = mace_mp(model="medium-mpa-0", device="cpu", default_dtype="float64", dispersion=False)

def sg(s, p):
    try: return SpacegroupAnalyzer(s, symprec=p).get_space_group_symbol()
    except Exception: return "?"

def prim(s):  # 24-atom conventional end members -> 12-atom primitive
    return SpacegroupAnalyzer(s, symprec=1e-3).get_primitive_standard_structure() if len(s) > 12 else s

# 1. reference structures
def make_ref(name, src, frm, to, note):
    d = CAND/name
    if (d/"structure.cif").exists(): return
    d.mkdir(exist_ok=True)
    s = Structure.from_file(CAND/src/"structure.cif"); s.replace_species({frm: to}) if False else None
    s = Structure.from_file(CAND/src/"structure.cif")
    s.replace_species({frm: to})
    s.to(filename=str(d/"structure.cif"))
    json.dump({"label": "PROPOSED MODEL - pure element substitution, unrelaxed, not computed",
               "source": src, "substitution": f"{frm}->{to}", "note": note,
               "formula": s.composition.reduced_formula}, open(d/"provenance.json", "w"), indent=2)
make_ref("REF-NaBiBr4-on-NaBiCl4proto", "NBCB-0-parent", "Cl", "Br", "all-Br on NaBiCl4 prototype (lattice not rescaled)")
make_ref("REF-NaBiCl4-on-NaBiBr4proto", "NBCB-4-NaBiBr4", "Br", "Cl", "all-Cl on NaBiBr4 prototype (lattice not rescaled)")

# 2. relaxations
def relax(cid):
    out = RES/cid
    if (out/"result.json").exists(): return
    out.mkdir(exist_ok=True)
    s0 = prim(Structure.from_file(CAND/cid/"structure.cif"))
    a = AA.get_atoms(s0); a.calc = calc
    t = time.time()
    f = FrechetCellFilter(a)
    opt = FIRE(f, logfile=None)
    conv = opt.run(fmax=0.01, steps=500)
    s1 = AA.get_structure(a)
    n = len(a)
    disp = max(np.linalg.norm(s1.lattice.get_cartesian_coords(
        (s1.frac_coords - s0.frac_coords + 0.5) % 1 - 0.5), axis=1)) if False else None
    # displacement: minimum-image Cartesian, in the relaxed cell
    df = (s1.frac_coords - s0.frac_coords + 0.5) % 1 - 0.5
    disp = float(np.max(np.linalg.norm(s1.lattice.get_cartesian_coords(df), axis=1)))
    dm = s1.distance_matrix + np.eye(n) * 1e6
    r = dict(id=cid, formula=s1.composition.reduced_formula, n_atoms=n, converged=bool(conv), steps=int(opt.nsteps),
             E_per_atom=float(a.get_potential_energy()/n), vol_per_atom_in=s0.volume/n, vol_per_atom=s1.volume/n,
             vol_change_pct=100*(s1.volume/s0.volume-1), sg_0p01=sg(s1, 0.01), sg_0p1=sg(s1, 0.1),
             max_disp_A=disp, min_dist_A=float(dm.min()), fmax_final=float(np.abs(a.get_forces()).max()),
             wall_s=time.time()-t)
    s1.to(filename=str(out/"relaxed.cif"))
    json.dump(r, open(out/"result.json", "w"), indent=2)
    print(cid, f"E={r['E_per_atom']:.5f} dV={r['vol_change_pct']:.1f}% conv={conv} steps={r['steps']} {r['sg_0p01']}", flush=True)

ids = sorted(p.name for p in CAND.iterdir() if p.is_dir())
for cid in ids: relax(cid)
R = {i: json.load(open(RES/i/"result.json")) for i in ids}

# 3. reference energies
cl = {k: R[k] for k in ["NBCB-0-parent", "REF-NaBiCl4-on-NaBiBr4proto"]}
br = {k: R[k] for k in ["NBCB-4-NaBiBr4", "REF-NaBiBr4-on-NaBiCl4proto"]}
kc = min(cl, key=lambda k: cl[k]["E_per_atom"]); kb = min(br, key=lambda k: br[k]["E_per_atom"])
ECl, EBr = R[kc]["E_per_atom"], R[kb]["E_per_atom"]
rows = []
for i in ids:
    r = R[i]
    if i.startswith("NBCB-") and i.split("-")[1] in "1234" and "-" in i[5:] and i.count("-") == 2 and i[5].isdigit():
        x = int(i.split("-")[1])
    elif i == "NBCB-0-parent": x = 0
    elif i == "NBCB-4-NaBiBr4": x = 4
    else: x = 0 if i.startswith("REF-NaBiCl4") else 4
    r["x"] = x
    r["E_mix_meV"] = 1000*(r["E_per_atom"] - ((4-x)/4*ECl + x/4*EBr))
    r["flags"] = ";".join(f for f, c in [("short_dist<1.8", r["min_dist_A"] < 1.8), ("dV>25%", abs(r["vol_change_pct"]) > 25),
                                          ("not_converged", not r["converged"])] if c)
    json.dump(r, open(RES/i/"result.json", "w"), indent=2)

# 4. phonons
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from ase import Atoms
def phonon(cid):
    out = PH/f"{cid}.json"
    if out.exists(): return json.load(open(out))
    s = Structure.from_file(RES/cid/"relaxed.cif")
    u = PhonopyAtoms(symbols=[str(x.specie) for x in s], cell=s.lattice.matrix, scaled_positions=s.frac_coords)
    ph = Phonopy(u, supercell_matrix=np.diag([2, 2, 2]), primitive_matrix=None)
    ph.generate_displacements(distance=0.01)
    t = time.time(); forces = []; budget = 600; aborted = False
    for sc in ph.supercells_with_displacements:
        at = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True); at.calc = calc
        forces.append(at.get_forces())
        if time.time()-t > budget: aborted = True; break
    res = dict(id=cid, n_disp=len(ph.supercells_with_displacements), done=len(forces), wall_s=time.time()-t)
    if aborted: res["aborted"] = "exceeded 600 s"
    else:
        ph.forces = forces; ph.produce_force_constants()
        ph.run_mesh([2, 2, 2], is_gamma_center=True, with_eigenvectors=False)
        f2 = ph.get_mesh_dict()["frequencies"]
        ph.run_mesh([8, 8, 8], is_gamma_center=True, with_eigenvectors=False)
        f8 = ph.get_mesh_dict()["frequencies"]
        res.update(min_THz_commensurate=float(f2.min()), min_THz_8x8x8=float(f8.min()),
                   n_modes_below_m0p1_commensurate=int((f2 < -0.1).sum()), n_modes_below_m0p1_8x8x8=int((f8 < -0.1).sum()),
                   n_qpoints_8x8x8=int(f8.shape[0]))
    json.dump(res, open(out, "w"), indent=2); print("phonon", res, flush=True)
    return res

sel = [kc, kb]
for x in (1, 2, 3):
    ok = sorted([r for r in R.values() if r["x"] == x and r["id"].startswith("NBCB-") and r["converged"] and not r["flags"]],
                key=lambda r: r["E_mix_meV"])[:3]
    sel += [r["id"] for r in ok]
P = {i: phonon(i) for i in sel}

# 5. outputs
tot = time.time()-T0
relax_wall = sum(r["wall_s"] for r in R.values())
ph_wall = sum(p["wall_s"] for p in P.values())
def pm(i):
    p = P.get(i); 
    if not p: return ""
    return f'{p["min_THz_8x8x8"]:.2f}' if "min_THz_8x8x8" in p else "aborted"
srt = sorted(R.values(), key=lambda r: r["E_mix_meV"])
cols = ["id", "x", "E_mix_meV", "vol_change_pct", "sg_0p01", "sg_0p1", "converged", "steps", "min_dist_A", "max_disp_A", "flags", "phonon_min_THz_8x8x8"]
with open(D/"results"/"screen_cycle2.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(cols)
    for r in srt: w.writerow([r.get(c, "") if c != "phonon_min_THz_8x8x8" else pm(r["id"]) for c in cols])
L = ["# Cycle 2 screen: MACE-MPA-0 (medium-mpa-0, float64, CPU, no dispersion)", "",
     "All structures are PROPOSED MODELS. MACE-MPA-0 mixing energies of a few meV/atom are within model error; "
     "these numbers only rank candidates for DFT and do not establish stability.", "",
     "## End members", "", "| structure | E (eV/atom) | dV vs input | SG(0.01) | converged |", "|---|---|---|---|---|"]
for k in list(cl)+list(br):
    r = R[k]; L.append(f'| {k} | {r["E_per_atom"]:.5f} | {r["vol_change_pct"]:+.1f}% | {r["sg_0p01"]} | {r["converged"]} |')
L += ["", f"- NaBiCl4 lowest-energy prototype (E_Cl reference): **{kc}** ; other prototype is {1000*(max(c['E_per_atom'] for c in cl.values())-ECl):.1f} meV/atom higher.",
      f"- NaBiBr4 lowest-energy prototype (E_Br reference): **{kb}** ; other prototype is {1000*(max(c['E_per_atom'] for c in br.values())-EBr):.1f} meV/atom higher.",
      f"- Sanity: NBCB-0-parent volume change from PBE-D3 input: {R['NBCB-0-parent']['vol_change_pct']:+.1f}%.",
      "- Reference-structure Br/Cl substitutions keep the donor lattice (not rescaled), so their reported dV is large by construction.", "",
      "## All structures sorted by E_mix", "",
      "E_mix = E - [(4-x)/4 E_Cl + x/4 E_Br], meV/atom. Phonon column = min frequency (THz) on 8x8x8 mesh, blank if not computed.", "",
      "| ID | x | E_mix (meV/atom) | dV % | SG 0.01 | SG 0.1 | converged | phonon min (THz) | flags |", "|---|---|---|---|---|---|---|---|---|"]
for r in srt:
    L.append(f'| {r["id"]} | {r["x"]} | {r["E_mix_meV"]:+.2f} | {r["vol_change_pct"]:+.1f} | {r["sg_0p01"]} | {r["sg_0p1"]} | {r["converged"]} | {pm(r["id"])} | {r["flags"]} |')
L += ["", "## Phonons (2x2x2 supercell of relaxed primitive cell, 0.01 A displacements)", "",
      "| ID | n disp | wall s | min THz commensurate | min THz 8x8x8 | modes < -0.1 THz (commens.) | modes < -0.1 THz (8x8x8) |", "|---|---|---|---|---|---|---|"]
for i, p in P.items():
    L.append(f'| {i} | {p["n_disp"]} | {p["wall_s"]:.0f} | {p.get("min_THz_commensurate","-")} | {p.get("min_THz_8x8x8","-")} | {p.get("n_modes_below_m0p1_commensurate","-")} | {p.get("n_modes_below_m0p1_8x8x8","-")} |'
             if "aborted" not in p else f'| {i} | {p["n_disp"]} | {p["wall_s"]:.0f} | aborted ({p["aborted"]}) | | | |')
n0 = sum(r["E_mix_meV"] <= 0 for r in R.values() if r["id"].startswith("NBCB-") and 1 <= r["x"] <= 3 and r["id"][5].isdigit() and r["id"].count("-") == 2)
n5 = sum(r["E_mix_meV"] <= 5 for r in R.values() if r["id"].startswith("NBCB-") and 1 <= r["x"] <= 3 and r["id"].count("-") == 2)
L += ["", f"Orderings (x=1-3, 42 total) with E_mix <= 0: {n0}; <= +5 meV/atom: {n5}.", "",
      "## Wall time", "", f"- Total script wall time this run: {tot/60:.1f} min (relaxations summed: {relax_wall/60:.1f} min; phonons: {ph_wall/60:.1f} min). Resumed results are not re-timed.", "",
      "## Caveats", "",
      "- MACE-MPA-0 is a PBE-like universal potential trained on MPtrj/Alexandria; no D3 dispersion here, so volumes tend to expand versus the PBE-D3 parent. Mixing energies of a few meV/atom are within model error.",
      "- The 24-atom conventional end-member cells were reduced to 12-atom primitive cells before relaxation (energy per atom is unchanged).",
      "- No symmetry constraints; FIRE can stall in saddle points of symmetric starts, and P1 orderings may be local minima only. Space groups are reported after relaxation at symprec 0.01 and 0.1.",
      "- Energies are 0 K total energies; no vibrational or configurational entropy. Phonons with 2x2x2 supercells sample q finely only on a Fourier interpolation; small negative values (> -0.1 THz) near Gamma can be numerical noise.",
      "- Nothing here is a stability or discovery claim; it is a ranking for DFT follow-up."]
open(D/"results"/"screen_cycle2.md", "w").write("\n".join(L)+"\n")
print("DONE", tot)
