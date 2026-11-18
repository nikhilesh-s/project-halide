# Cycle 2 screen: MACE-MPA-0 (medium-mpa-0, float64, CPU, no dispersion)

All structures are PROPOSED MODELS. MACE-MPA-0 mixing energies of a few meV/atom are within model error; these numbers only rank candidates for DFT and do not establish stability.

## End members

| structure | E (eV/atom) | dV vs input | SG(0.01) | converged |
|---|---|---|---|---|
| NBCB-0-parent | -3.37487 | +11.8% | C2/c | True |
| REF-NaBiCl4-on-NaBiBr4proto | -3.35323 | -13.5% | C2/c | True |
| NBCB-4-NaBiBr4 | -2.98099 | +0.9% | C2/c | True |
| REF-NaBiBr4-on-NaBiCl4proto | -2.99642 | +27.5% | C2/c | True |

- NaBiCl4 lowest-energy prototype (E_Cl reference): **NBCB-0-parent** ; other prototype is 21.6 meV/atom higher.
- NaBiBr4 lowest-energy prototype (E_Br reference): **REF-NaBiBr4-on-NaBiCl4proto** ; other prototype is 15.4 meV/atom higher.
- Sanity: NBCB-0-parent volume change from PBE-D3 input: +11.8%.
- Reference-structure Br/Cl substitutions keep the donor lattice (not rescaled), so their reported dV is large by construction.

## All structures sorted by E_mix

E_mix = E - [(4-x)/4 E_Cl + x/4 E_Br], meV/atom. Phonon column = min frequency (THz) on 8x8x8 mesh, blank if not computed.

| ID | x | E_mix (meV/atom) | dV % | SG 0.01 | SG 0.1 | converged | phonon min (THz) | flags |
|---|---|---|---|---|---|---|---|---|
| NBCB-0-parent | 0 | +0.00 | +11.8 | C2/c | C2/c | True | -0.71 |  |
| REF-NaBiBr4-on-NaBiCl4proto | 4 | +0.00 | +27.5 | C2/c | C2/c | True | -0.70 | dV>25% |
| NBCB-3-08 | 3 | +0.32 | +22.7 | P-1 | P-1 | True | -0.64 |  |
| NBCB-2-06 | 2 | +0.49 | +19.5 | P-1 | P-1 | True | -0.63 |  |
| NBCB-3-05 | 3 | +1.01 | +23.9 | P1 | P1 | True | -0.55 |  |
| NBCB-1-01 | 1 | +1.23 | +16.3 | P-1 | P-1 | True | -0.71 |  |
| NBCB-2-08 | 2 | +1.41 | +19.4 | P1 | P1 | True | -0.52 |  |
| NBCB-3-10 | 3 | +1.44 | +21.5 | Cc | Cc | True | -0.61 |  |
| NBCB-3-07 | 3 | +1.45 | +23.3 | P1 | P1 | True |  |  |
| NBCB-1-10 | 1 | +1.46 | +14.7 | Cc | Cc | True | -0.57 |  |
| NBCB-1-08 | 1 | +1.51 | +14.5 | P-1 | P-1 | True | -0.83 |  |
| NBCB-1-05 | 1 | +1.51 | +15.7 | P1 | P1 | True |  |  |
| NBCB-2-19 | 2 | +1.63 | +18.7 | P1 | P1 | True | -0.69 |  |
| NBCB-2-07 | 2 | +1.71 | +19.8 | P1 | P1 | True |  |  |
| NBCB-2-22 | 2 | +1.73 | +17.1 | C2/c | C2/c | True |  |  |
| NBCB-3-01 | 3 | +1.77 | +24.3 | P-1 | P-1 | True |  |  |
| NBCB-3-04 | 3 | +1.82 | +23.4 | P1 | P1 | True |  |  |
| NBCB-1-07 | 1 | +1.85 | +15.4 | P1 | P1 | True |  |  |
| NBCB-1-09 | 1 | +1.96 | +15.0 | C2 | C2 | True |  |  |
| NBCB-2-09 | 2 | +2.00 | +19.3 | P-1 | P-1 | True |  |  |
| NBCB-1-04 | 1 | +2.05 | +15.3 | P1 | P1 | True |  |  |
| NBCB-3-09 | 3 | +2.08 | +22.5 | C2 | C2 | True |  |  |
| NBCB-3-03 | 3 | +2.09 | +24.3 | Cc | Cc | True |  |  |
| NBCB-2-21 | 2 | +2.13 | +17.9 | P1 | P1 | True |  |  |
| NBCB-3-06 | 3 | +2.15 | +23.6 | P1 | P1 | True |  |  |
| NBCB-2-13 | 2 | +2.17 | +19.6 | C2 | C2 | True |  |  |
| NBCB-2-18 | 2 | +2.18 | +18.2 | P1 | P1 | True |  |  |
| NBCB-2-03 | 2 | +2.22 | +19.8 | P1 | P1 | True |  |  |
| NBCB-3-02 | 3 | +2.24 | +25.1 | C2 | C2 | True |  | dV>25% |
| NBCB-2-05 | 2 | +2.35 | +20.3 | P1 | P1 | True |  |  |
| NBCB-1-06 | 1 | +2.39 | +15.5 | P1 | P1 | True |  |  |
| NBCB-2-14 | 2 | +2.40 | +19.0 | P1 | P1 | True |  |  |
| NBCB-2-17 | 2 | +2.47 | +19.3 | Cc | Cc | True |  |  |
| NBCB-1-03 | 1 | +2.47 | +15.8 | Cc | Cc | True |  |  |
| NBCB-2-16 | 2 | +2.56 | +19.0 | Cc | Cc | True |  |  |
| NBCB-2-20 | 2 | +2.63 | +18.0 | P1 | P1 | True |  |  |
| NBCB-2-02 | 2 | +2.83 | +21.0 | P1 | P1 | True |  |  |
| NBCB-2-15 | 2 | +2.84 | +19.5 | P1 | P1 | True |  |  |
| NBCB-2-10 | 2 | +3.02 | +19.4 | P1 | P1 | True |  |  |
| NBCB-2-01 | 2 | +3.05 | +21.0 | C2/c | C2/c | True |  |  |
| NBCB-2-12 | 2 | +3.10 | +19.8 | P1 | P1 | True |  |  |
| NBCB-1-02 | 1 | +3.11 | +17.2 | C2 | C2 | True |  |  |
| NBCB-2-04 | 2 | +3.25 | +20.8 | P1 | P1 | True |  |  |
| NBCB-2-11 | 2 | +4.89 | +19.5 | C2 | C2 | True |  |  |
| NBCB-4-NaBiBr4 | 4 | +15.43 | +0.9 | C2/c | C2/c | True |  |  |
| REF-NaBiCl4-on-NaBiBr4proto | 0 | +21.64 | -13.5 | C2/c | C2/c | True |  |  |

## Phonons (2x2x2 supercell of relaxed primitive cell, 0.01 A displacements)

| ID | n disp | wall s | min THz commensurate | min THz 8x8x8 | modes < -0.1 THz (commens.) | modes < -0.1 THz (8x8x8) |
|---|---|---|---|---|---|---|
| NBCB-0-parent | 18 | 5 | -0.00019525612893009266 | -0.709043884268494 | 0 | 83 |
| REF-NaBiBr4-on-NaBiCl4proto | 18 | 5 | -0.17823474578523768 | -0.6955136594113309 | 1 | 81 |
| NBCB-1-01 | 36 | 10 | -0.00010557165677542251 | -0.7069710990676453 | 0 | 122 |
| NBCB-1-10 | 36 | 10 | -0.0003154503413047787 | -0.5719406663030702 | 0 | 66 |
| NBCB-1-08 | 36 | 10 | -8.939727043435302e-05 | -0.8329175089729816 | 0 | 156 |
| NBCB-2-06 | 36 | 9 | -6.901056428098711e-05 | -0.6283779923190083 | 0 | 140 |
| NBCB-2-08 | 72 | 19 | -8.376337751249927e-05 | -0.5224366567082706 | 0 | 138 |
| NBCB-2-19 | 72 | 20 | -6.995146447177884e-05 | -0.6869525226004176 | 0 | 164 |
| NBCB-3-08 | 36 | 11 | -6.368961963895922e-05 | -0.6445679901387058 | 0 | 150 |
| NBCB-3-05 | 72 | 20 | -7.943150485578151e-05 | -0.5524412444142979 | 0 | 135 |
| NBCB-3-10 | 36 | 10 | -0.0002673158956416095 | -0.612936368201842 | 0 | 87 |

Orderings (x=1-3, 42 total) with E_mix <= 0: 0; <= +5 meV/atom: 42.

## Wall time

- Total script wall time this run: 0.0 min (relaxations summed: 2.2 min; phonons: 2.1 min). Resumed results are not re-timed.

## Caveats

- MACE-MPA-0 is a PBE-like universal potential trained on MPtrj/Alexandria; no D3 dispersion here, so volumes tend to expand versus the PBE-D3 parent. Mixing energies of a few meV/atom are within model error.
- The 24-atom conventional end-member cells were reduced to 12-atom primitive cells before relaxation (energy per atom is unchanged).
- No symmetry constraints; FIRE can stall in saddle points of symmetric starts, and P1 orderings may be local minima only. Space groups are reported after relaxation at symprec 0.01 and 0.1.
- Energies are 0 K total energies; no vibrational or configurational entropy. Phonons with 2x2x2 supercells sample q finely only on a Fourier interpolation; small negative values (> -0.1 THz) near Gamma can be numerical noise.
- Nothing here is a stability or discovery claim; it is a ranking for DFT follow-up.
- Phonon interpretation: at the commensurate q-points (Gamma and zone-boundary points of the 2x2x2 supercell) the minima are ~0 THz (acoustic noise) for every computed structure except REF-NaBiBr4-on-NaBiCl4proto (-0.18 THz, 1 mode). The 8x8x8 minima (-0.5 to -0.8 THz, 60-160 modes below -0.1 THz) come from Fourier interpolation of 2x2x2 force constants and are not reliable evidence for or against dynamical stability; a larger supercell or DFT is needed.
- The lowest-energy NaBiBr4 reference is the NaBiCl4-prototype substitution (a proposed model, not a known structure; dV +27.5% flagged), so E_Br and every E_mix depends on it; against the real NaBiBr4 prototype, E_mix would shift by -(x/4)*15.4 meV/atom (i.e. lower).
- Timing: the first full run took ~4 min total (relaxations ~2 min, phonons 5-20 s per job, all far under the 10 min cap); the totals above are from the final resumed run.
