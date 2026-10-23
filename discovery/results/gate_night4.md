# Gate night 4: which remaining cycle-3 hits get DFT

Every TP-* structure below is a **PROPOSED MODEL** from a prototype transplant, ranked by MACE-MPA-0 + D3(BJ). Nothing here is a stability or discovery claim. No DFT has been run on any structure in this file.

Method (identical to gate 3): best hit per chemistry from `transplant_cycle3.csv` (lowest decomposition energy among hits with dE <= -5 meV/atom, primitive cell <= 12 atoms). The hit, the three lowest Alexandria AMX4 entries by hull distance (<= 16 atoms) plus the scan's best baseline, and every Alexandria AX / MX3 / ternary entry with hull distance <= 5 meV/atom and <= 20 atoms were rattled 0.03 Å (seed 3) and re-relaxed: MACE-MPA-0 float64 + D3, FIRE + FrechetCellFilter, fmax 0.01, 3 single-thread workers. 78 relaxations, 27 min, all converged. Keep rule: dE <= -5 and decomposition vs AX + MX3 <= -3 meV/atom after rattle. Scripts: `discovery/scripts/gate4/{run,analyze,make_night4}.py` (`run.py` is gate3's with the chemistry list changed); raw: `discovery/results/gate4/`.

Chemistries checked: the seven with a cycle-3 hit passing both cuts that are not in any DFT queue (night2, night3, air_day). Na-In, Na-Sc, NaBiI4, NaSbI4 hits have positive decomposition in the scan and were not re-checked. The Li iodides and NaInI4/NaScI4/NaYI4 were never scanned in cycle 3 (dropped for budget), so there are no iodide hits other than NaGaI4.

## ML gate result

| Hit | Atoms | dE vs best existing, before -> after (meV/atom) | Decomp vs AX + MX3, before -> after | Best existing after rattle (others, meV above it) | ML gate | DFT queue |
|---|---|---|---|---|---|---|
| TP-LiInCl4-01 | 12 | -7.0 -> -7.0 | -10.0 -> -8.0 | agm003740637 (+50.5, +54.2) | pass | yes (01-04) |
| TP-LiScBr4-01 | 12 | -6.8 -> -6.8 | -10.9 -> -10.8 | agm005739333 (+43.3, +122.1) | pass | yes (05-09) |
| TP-LiYBr4-02 | 12 | -7.0 -> -7.0 | -7.9 -> -7.2 | agm003742437 (+67.6) | pass | yes (10-12) |
| TP-NaGaI4-02 | 12 | -18.1 -> -18.2 | -5.0 -> -5.0 | agm003736187 (+17.9 agm005755304, +31.3) | pass | yes, last (13-17) |
| TP-LiInBr4-01 | 12 | -7.5 -> -7.4 | -11.4 -> -11.4 | agm003740511 (+37.9, +39.9) | pass | **no**: experimental structure known (disordered Fd-3m spinel, the literature notes) |
| TP-LiAlBr4-01 | 12 | -16.3 -> -16.2 | -47.9 -> -48.0 | agm003740512 (+1.7, +7.7) | pass (gate 3) | **no**: experimental structure known (P2_1/c, Z = 4) |
| TP-LiScCl4-02 | 12 | -11.0 -> -11.0 | -0.8 -> -0.8 | agm005742973 (+23.1, +52.7) | **fail** (decomposition > -3) | no |

Four kept for DFT, below the 5-7 target: after the scout's known-structure list removes LiAlBr4 and LiInBr4 there is nothing else in the cycle-3 hit list that passes. More kept chemistries need a new scan (Li iodides, Na-M-I), not a looser gate.

For every chemistry the lowest MACE+D3 decomposition route among all relaxed competitors is AX + MX3 (the listed ternaries LiSc3Br8, LiSc3Cl8, NaGaI3 do not lower the hull at the 1:1:4 composition), so the full-hull number equals the binary one.

Caveats, in order of importance:
1. **The margins are small and D3-dependent.** In PBE without D3, Alexandria places the existing LiInCl4 / LiScBr4 / LiYBr4 entries 13.1 / 18.5 / 26.9 meV/atom above its hull; the hits are only 7 below those entries. The negative decomposition energies come from MACE+D3 and are within its error (RMS 3-11 meV/atom in the systems tested so far). Expect some of these to fail DFT.
2. **Missing competitors.** Li3InCl6, Li3YBr6 and Li3ScCl6-type phases are real, synthesised compounds, but their Alexandria entries are 26-28 meV/atom above the hull (Li3InCl6 agm003569151 +28.3, Li3YBr6 agm003224330 +26.0), outside the 5 meV rule, so they are not in the queue. the literature notes says any LiInCl4 claim must beat Li3InCl6 + InCl3; a DFT pass on jobs 01-04 is therefore G2 plus binaries only, not G3.
3. **Not scouted:** LiScBr4, LiYBr4, NaGaI4 were outside the literature notes's scope. NaGaI4 has three Alexandria entries within 4 meV of the hull and its Cl/Br analogues have known structures, so it is the likeliest to be a known compound; it is queued last.
4. NaGaI4: the Alexandria hull-0 entry (agm005755304, 6 atoms) is 17.9 meV/atom above agm003736187 in MACE+D3; both are in the queue.

## Method validation: does the search recover known experimental structures? (ML only, no DFT inputs written)

Experimental space groups from the literature notes. "Our lowest" = lowest MACE+D3 energy among all Alexandria entries (<= 24 atoms) and all transplants for that composition in the cycle-3 scan.

| Compound | Experimental | Listed in Alexandria with that space group? | Our lowest MACE+D3 structure | Recovered? |
|---|---|---|---|---|
| NaAlBr4 | Pnma, Z = 4 | agm003232517 (SG 62, 24 atoms, hull 0) | agm003232517, Pnma; best transplant -0.01 | **yes** |
| NaGaCl4 | P2_12_12_1, Z = 4 | agm003232627 (SG 19, 24 atoms, hull 0) | agm003232627, P2_12_12_1; best transplant +4.4 | **yes** |
| LiGaCl4 | P2_1/c, Z = 4 | agm003232631 (SG 14, 24 atoms, PBE hull +4.9; PBE prefers SG 82) | agm003232631, P2_1/c; the lowest transplant (LiGaI4-type, -7.6) is a loose match to the same entry | **yes** (and D3 fixes the PBE ordering) |
| LiAlBr4 | P2_1/c, Z = 4 | **no** SG 14 entry (listed: SG 148 36 atoms hull 0, SG 82, 112, 20) | TP-LiAlBr4-01, P2_1/m, 12 atoms, -16.3; the LiAlCl4-type P2_1/c transplant (24 atoms, from agm003231472) is +3.2 vs the best listed, about 19 above TP-01; lowest P2_1/c transplant is KGaCl4-type at -10.7 | **no** |

Reading: three of four recovered by space group and cell size (no atom-by-atom match against experimental coordinates was done; none are in the repo). LiAlBr4 is a miss: MACE+D3 ranks a P2_1/m model about 19 meV/atom below the experimental structure type. Either the potential is wrong for Li-Al-Br by that much, or the model is a 0 K polymorph; this cannot be settled at ML level, and it is the same size as the dE that carries the four kept hits. Also not relaxed: the 36-atom Alexandria LiAlBr4 entry agm081834437 (hull 0), above the size cap in gate 3 and in the scan. It was assumed, not checked, that agm003231472 is the experimental LiAlCl4 structure.

## Cutoffs

Rule: ecutwfc = max(60, highest "suggested minimum cutoff for wavefunctions" among the reaction's elements, rounded up to 60/80/100 Ry), ecutrho = 8x, one cutoff for every phase of a reaction. Header values (Ry): Li 103, Na 66, Ga 60, Sc 50, Br 49, Cl 45, Bi 45, In 41, Y 38, I 36, Al 29, Sb 27.

| Reaction | Driver | ecutwfc / ecutrho |
|---|---|---|
| Li-In-Cl, Li-Sc-Br, Li-Y-Br | Li 103 (capped at the top tier) | 100 / 800 |
| Na-Ga-I | Na 66 | 80 / 640 |

Two consequences. (a) This is stricter than earlier nights (Li at 80, Na at 60), so the earlier LiCl / LiBr (80 Ry) and NaCl (60 Ry) results cannot be reused and LiCl, LiBr are recomputed at 100 Ry. (b) Li at 100 Ry is slightly below the 103 Ry suggestion. Running the Li jobs at 80 Ry instead would match nights 2-3, cut about 45% of the time and allow reuse of LiBr from night 3; that is a one-line change in `make_night4.py` (`SUG`/tiers) if preferred. ecutrho 8x covers every charge-density suggestion (largest: Bi 455, Li 410, Sc 402).

## DFT queue: `overnight/night4/jobs/` (17 jobs, about 22.6 h, i.e. about three nights)

Settings as nights 1-3 (`make_inputs.write`, vc-relax, PBE-D3, k-spacing 0.22). Starting geometries: the rattled, re-relaxed cells in `discovery/results/gate4/relax/`. Estimates scale night2 (12 atoms at 60 Ry 31-73 min) by 2 at 80 Ry and 3.5 at 100 Ry; they are rough, and the In/Ga (13 valence electrons) and Sc/Y (semicore) cells have not been timed. Atom counts are primitive cells as written.

| File (suffix `_vcrelax.in`) | Atoms | Ecut (Ry) | Est. (min) |
|---|---|---|---|
| 01_LiInCl4_TP01_hit | 12 | 100 | 150 |
| 02_LiInCl4_agm003740637_best_existing | 6 | 100 | 50 |
| 03_LiCl_wz_agm002157414_100Ry_ref | 4 | 100 | 15 |
| 04_InCl3_agm002142774_ref | 8 | 100 | 70 |
| 05_LiScBr4_TP01_hit | 12 | 100 | 150 |
| 06_LiScBr4_agm005739333_best_existing | 6 | 100 | 50 |
| 07_LiBr_agm002091979_100Ry_ref | 4 | 100 | 15 |
| 08_ScBr3_agm002184240_ref | 8 | 100 | 70 |
| 09_LiSc3Br8_agm006022042_ternary (Alexandria hull 0) | 12 | 100 | 150 |
| 10_LiYBr4_TP02_hit | 12 | 100 | 150 |
| 11_LiYBr4_agm003742437_best_existing | 6 | 100 | 50 |
| 12_YBr3_agm003297005_ref | 8 | 100 | 70 |
| 13_NaGaI4_TP02_hit | 12 | 80 | 90 |
| 14_NaGaI4_agm003736187_best_existing | 12 | 80 | 90 |
| 15_NaGaI4_agm005755304_alex_hull0 | 6 | 80 | 30 |
| 16_NaI_agm003157588_ref | 2 | 80 | 5 |
| 17_GaI3_agm003220383_ref | 16 | 80 | 150 |

Per chemistry: LiInCl4 4.75 h, LiScBr4 7.25 h, LiYBr4 4.5 h (uses LiBr from job 07), NaGaI4 6.1 h. Overflow, `overnight/night4/overflow/` (not auto-run): 18_NaGaI3_agm005838860_ternary (20 atoms, 80 Ry, about 200 min; Alexandria hull 0 but Ga(II), not on the NaI-GaI3 line). The queue is one over the 16-job cap because both NaGaI4 baselines are kept. Only the three lowest listed AMX4 entries were ML-checked per chemistry and one or two are in DFT; the others are 18-120 meV/atom higher in MACE+D3.

## Pseudopotentials added (PSlibrary 1.0.0 PBE PAW, scalar-relativistic, downloaded to `overnight/pseudo/`, appended to `overnight/pseudo_md5.txt`)

| File | MD5 | Suggested ecutwfc / ecutrho (Ry) | Valence e | Used by night 4 |
|---|---|---|---|---|
| In.pbe-dn-kjpaw_psl.1.0.0.UPF | 235b15e6a7fdfd96486e0ecdcb5b12c6 | 41 / 212 | 13 | yes |
| Sc.pbe-spn-kjpaw_psl.1.0.0.UPF | d5ac433e587a0148d610ca85a1be5559 | 50 / 402 | 11 | yes |
| Y.pbe-spn-kjpaw_psl.1.0.0.UPF | d52e8d8ffa9ef5a9f13bb13aaba2e897 | 38 / 257 | 11 | yes |
| Ga.pbe-dn-kjpaw_psl.1.0.0.UPF | 50ce49cbf7e0a40ad8623bf41840a268 | 60 / 244 | 13 | yes |
| I.pbe-n-kjpaw_psl.1.0.0.UPF | cd3d296c3308a581e1559a0c8a554cb5 | 36 / 154 | 7 | yes |
| Al.pbe-n-kjpaw_psl.1.0.0.UPF | 502b6b27b9c6270cab5dd1084b96fa6d | 29 / 143 | 3 | no (LiAlBr4 dropped) |

`setup.sh` must be re-run once on wifi before night 4 (about 11 MB). Dry start: `mpirun -np 2 pw.x` for 40 s each on jobs 04 (In), 08 (Sc), 12 (Y), 17 (Ga, I) with placeholders filled in a temp copy: all pseudopotentials read with the MD5s above, cutoffs and k-points printed, no input error. None reached the first SCF iteration in 40 s, so SCF convergence with these files is untested. Temp files removed.

Candidate files: `discovery/candidates/<ID>/gate4_check.json` and `relaxed_mace_d3_rattled_gate4.cif` for the four kept hits.
