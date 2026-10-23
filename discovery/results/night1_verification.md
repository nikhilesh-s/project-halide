# Night 1 verification (read-only; script: discovery/scripts/verify_night1.py)

Energies = last "!" line (final SCF, fresh basis). Ry->meV: 13605.693123. 12 atoms (2 f.u.) for all NaBiCl4.

## 1. Output checks
| job | JOB DONE | bfgs converged | "NOT achieved"/Error | bfgs steps | last Total force (Ry/bohr) | final P (kbar) | last "!" (Ry) | Final enthalpy (Ry) | E(!) - H |
|---|---|---|---|---|---|---|---|---|---|
| 00 control (scf) | yes | n/a (scf) | none | - | 0.002087 | -1.13 | -2667.12492286 | - | - |
| ref NaBiCl4_d3 | yes | yes (34 steps) | none | 34 | 0.002132 | -1.12 | -2667.12492286 | -2667.13299844 | +8.08 mRy = +9.16 meV/atom |
| 01 agm088474086 | yes | yes (16 steps, 17 scf) | none | 16 | 0.001226 (relax) | -0.10 | -2667.14146589 | -2667.14153071 | +0.07 mRy = +0.07 meV/atom |
| 02 agm084661913 | yes | yes (11 steps, 12 scf) | none | 11 | 0.001645 (relax) | +0.23 | -2667.13805369 | -2667.13769469 | -0.36 mRy = -0.41 meV/atom |

- Control 00 reproduces -2667.12492286 Ry exactly (identical to reference "!"). PASS.
- NOTE: the reference's final-enthalpy (variable-basis, last bfgs step) is 9.2 meV/atom BELOW its fresh-basis "!" energy. The reference started at 2800.6 a.u.^3 and shrank to 2468.1 a.u.^3 (-12%), so its relaxed-basis enthalpy suffers Pulay-stress artefact. Fresh-basis "!" is the right quantity (same fixed cutoff for all), but if one had compared Final enthalpies the gap would be only 9.7 (01) / 5.3 (02) meV/atom. For 01/02 the two energies agree to <0.5 meV/atom, so only the reference is affected. 60 Ry/480 Ry is thus not fully converged for volume changes; absolute energies are cutoff-dependent but the fresh-basis comparison at fixed cutoff is consistent.
- Residual P at the fresh-basis final SCF: ref -1.12 kbar, so the reference is close to equilibrium at this cutoff; the 9 meV/atom basis shift is not an unrelaxed-geometry artefact.

## 2. Settings vs reference (all identical unless noted)
ecutwfc 60 / ecutrho 480 Ry; vdw_corr dft-d3 (3-body printed in all); Gaussian smearing degauss 0.005 Ry; conv_thr 1e-8, forc_conv_thr 1e-3, etot_conv_thr 1e-5, press_conv_thr 0.5; QE v.7.5; 12 atoms; 104 electrons.
Pseudos: Na.pbe-spn-kjpaw_psl.1.0.0.UPF (md5 05b4deb74e3932c23ab016979dc8b5d4), Cl.pbe-n-kjpaw_psl.1.0.0.UPF (530da96d25dcb7781b0586ca00827062), Bi.pbe-dn-kjpaw_psl.1.0.0.UPF (266ded784ae3bc462508ed03f0f37963) - same names and MD5 in 00, 01, 02, ref, NaCl, BiCl3 (NaCl lacks Bi; BiCl3 lacks Na). PASS.
K-meshes (Gamma-centred-equivalent 0 0 0 shift, automatic): ref and 00 = 4x5x4 (42 irreducible k); 01 = 5x4x5 (39 k); 02 = 5x5x5 (63 k, different cell/setting); NaCl 8x8x8; BiCl3 4x4x3.
K-mesh caveat (UNTESTED): all cells are 6-10 A per edge, wide-gap ionic halides, Gaussian smearing 0.005 Ry on an insulator (effectively no smearing). Typical k-convergence error for such systems at 4-5 divisions is of order 1 meV/atom or below, so an error of 15-19 meV/atom from k-sampling alone is implausible, but I did not test it; densities (k*a) are comparable between jobs (ref ~4x6.9..., 01 5x6.2...), so differential error should be a few meV/atom at most. A cheap SCF re-run of 01 at the 4x5x4-equivalent density and of ref at 5x5x5 would settle it.

## 3. Energy differences (recomputed)
- 01 - ref = (-2667.14146589 - -2667.12492286) = -0.01654303 Ry = -225.1 meV/cell = **-18.76 meV/atom** (matches claim 18.8)
- 02 - ref = -0.01313083 Ry = -178.6 meV/cell = **-14.89 meV/atom** (matches claim 14.9)
- 02 - 01 = +3.87 meV/atom (01 lowest)
- Using Final enthalpies instead: 01 -9.67, 02 -5.33 meV/atom (see Pulay note above).

## 4. Structures (final blocks, crystal coordinates; Angstrom cell)
| | ref (005743614 relaxed) | 01 (088474086) | 02 (084661913) |
|---|---|---|---|
| space group symprec 0.01 / 0.1 | C2/c (15) / C2/c | P2/c (13) / P2/c | Fdd2 (43) / Fdd2 |
| volume (A^3), per atom | 365.74, 30.48 | 336.12, 28.01 | 347.81, 28.98 |
| Bi-Cl CN (cut 3.3 A) | 6, 6 | 6, 6 | 6, 6 |
| Na-Cl CN (cut 3.3 A) | 6, 6 | 6, 6 | 6, 6 |
| shortest Bi-Cl (A) | 2.617 | 2.607 | 2.609 |
StructureMatcher (ltol .2/stol .3/angle 5 and loose .4/.6/10): ref-01 False, ref-02 False, 01-02 False under both. Three distinct structures; neither 01 nor 02 relaxed into the reference or each other (different space groups, volumes 8.1% and 5.2% smaller than ref... 01 is 8.1% smaller, 02 4.9% smaller).
Note: all three have 6-coordinate Bi and Na, so differences are in packing/connectivity, not coordination; ref has the largest volume/atom by 9%, consistent with the energy ordering E(V) at fixed cutoff.

## 5. Decomposition vs NaCl + BiCl3 (same cutoffs, pseudos, D3 confirmed)
NaCl_d3.out: 2 atoms = 1 f.u., last "!" -188.87219971 Ry. BiCl3_d3.out: 16 atoms = 4 f.u., last "!" -4578.74871153 Ry (-1144.68717788 Ry/f.u.). Reference mix per NaBiCl4 f.u. = -1333.55937759 Ry. Per-atom (6 atoms/f.u.), negative = stable:
- ref: **-6.99 meV/atom** (earlier report -7 to -9: reproduced)
- 01: **-25.75 meV/atom**
- 02: **-21.88 meV/atom**
Caveat: the BiCl3 fresh-basis "!" is 10.2 mRy (~8.7 meV/atom) above its own final enthalpy, same Pulay effect, so decomposition values depend on the basis-set treatment; all use consistent fresh-basis "!" here. NaCl and BiCl3 used different k meshes (8x8x8, 4x4x3), no k-test.

## Verdict
1. Numbers are internally consistent: control reproduces, all jobs converged, identical functional/cutoff/pseudos/D3/QE version, differences recomputed (-18.76, -14.89 meV/atom).
2. Structures 01 and 02 are genuinely distinct from the reference and from each other (P2/c, Fdd2 vs C2/c), no relaxation collapse.
3. Claim "088474086 lowest NaBiCl4 found so far at PBE-D3, ~19 meV/atom below" is SUPPORTED for fresh-basis energies at 60/480 Ry, within the polymorphs examined (only 3 tried).
4. Caveats: k-meshes differ between jobs (untested, likely a few meV/atom or less); reference suffers a 9 meV/atom Pulay shift (ref started 12% too large in volume; on enthalpies the gap is ~9.7 / ~5.3); single functional, no phonons/ZPE/finite-T, cutoff convergence for volume change not checked.
5. Before using it in a paper: rerun ref (or SCF at final geometry) at a higher cutoff (e.g. 80/640) and at a common k-mesh for all three; the ordering likely holds but the ~19 magnitude may shrink.
