# Gate night 3: which cycle-3 hits get DFT

Every TP-* structure below is a **PROPOSED MODEL** from a prototype transplant, ranked by MACE-MPA-0 + D3(BJ). Nothing here is a stability or discovery claim; the only DFT facts are those already in `overnight/night2` (TP-NaBiCl4-01 and TP-NaSbCl4-01 survived; TP-LiBiCl4-01 still running).

Method: best hit per chemistry from `transplant_cycle3.csv` (lowest decomposition energy among hits with dE <= -5 meV/atom); chemistries with best decomposition > -3 dropped; eight chemistries checked. Each hit and every Alexandria AMX4 entry of that chemistry with <= 16 atoms (the three lowest by Alexandria hull distance plus the scan's best baseline) was rattled 0.03 Å (seed 3) and re-relaxed, MACE-MPA-0 float64 + D3, FIRE + FrechetCellFilter, fmax 0.01, 3 CPU threads. Binaries: lowest-energy AX and MX3 among Alexandria entries with hull <= 5 meV/atom (relaxed the same way). Decomposition = E(AMX4) - [2E(AX) + 4E(MX3)]/6 per atom (negative = stable against the binaries). NaSbCl4, NaBiCl4, LiBiCl4 reuse the identical-setting skeptic3 relaxations. Scripts: `discovery/scripts/gate3/{run,analyze,make_night3}.py`; raw: `discovery/results/gate3/`.

Substitutions to the brief: eight chemistries are the top-8 by decomposition energy except LiInBr4 (needs In pseudopotential) was swapped for NaSbBr4 (DFT-feasible). LiSbCl4: the best hit (TP-LiSbCl4-01) has a 24-atom primitive cell, over the 16-atom ML cap, so the best hit with <= 12 atoms (TP-LiSbCl4-03) is used. LiBr reference uses rock salt (2 atoms) instead of the 4-atom entry that is 0.3 meV/atom lower in MACE.

## Kept hits (dE <= -5 and decomposition <= -3 after rattle)

| Hit | Atoms | dE vs best existing, before -> after (meV/atom) | Decomp vs AX + MX3, before -> after | Best existing after rattle | Status |
|---|---|---|---|---|---|
| TP-NaSbCl4-01 | 12 | -38.0 -> -9.9 | -18.7 -> -18.1 (-24.6 vs 1/3 Na3SbCl6 + 2/3 SbCl3) | agm005743441 | DFT survived; needs Na3SbCl6 |
| TP-LiBiCl4-01 | 12 | -33.6 -> -27.8 | -13.6 -> -6.9 (D3 only) | agm003736573 | DFT running (night2) |
| TP-LiSbBr4-01 | 12 | -23.7 -> -23.8 | -19.7 -> -19.8 | agm005739369 | queued |
| TP-LiBiBr4-01 | 12 | -23.4 -> -23.4 | -11.1 -> -11.1 | agm005739442 | queued |
| TP-NaSbBr4-01 | 12 | -22.0 -> -22.0 | -5.1 -> -5.2 | agm003733885 | overflow (next night) |
| TP-LiSbCl4-03 | 12 | -64.9 -> -46.3 | -6.4 -> -5.9 | agm005453208 (was agm005743014, now 18.7 higher) | overflow (next night) |
| TP-LiAlBr4-01 | 12 | -16.3 -> -16.2 | -47.9 -> -48.0 | agm003740512 | needs Al pseudopotential |

Caveats: ML trust in Sb systems is poor (MACE vs PBE RMS 11.4 meV/atom in Na-Sb-Cl, from skeptic3); the LiBiCl4 decomposition margin exists only with D3. The bromide numbers did not move because their baselines did not slide on rattle, which is reassuring but is not DFT.

## Rejected

| Chemistry | Hit | Reason |
|---|---|---|
| NaBiCl4 | TP-NaBiCl4-01 | After rattle dE = -2.1 (was -8.1 in the scan), above the -5 cut; baseline agm084661913 dropped 5.6. Already DFT-tested. |
| NaBiBr4 | TP-NaBiBr4-01 | Best decomposition -0.9 (> -3). |
| NaInCl4, NaInBr4, NaScBr4, NaBiI4, NaSbI4 | TP-*-01 | Decomposition > 0 (+1.5, +18.6, +12.2, +5.7, +6.2). |
| LiInBr4, LiScBr4, LiYBr4, LiInCl4, LiScCl4, NaGaI4 | TP-*-01 | Pass the decomposition cut but beyond the 8-chemistry cap and need new pseudopotentials; not rattle-checked. |

## DFT queue: `overnight/night3/jobs/` (8 jobs, about 5.9 h)

Same settings as night 1/2 (`make_inputs.write`, vc-relax, D3). Estimates from `night2/summary.json` (12 atoms at 60 Ry 31-73 min; the 12-atom Li job at 80 Ry took 71 min) with Br assumed slightly heavier.

| File (suffix `_vcrelax.in`) | Atoms | Ecut (Ry) | Est. |
|---|---|---|---|
| 01_Na3SbCl6_agm017656109 (highest priority; Alexandria hull +2.7 meV; Materials Project also lists it, not checked here) | 10 | 60 | 40 min |
| 02_LiSbBr4_TP01_hit | 12 | 80 | 90 min |
| 03_LiSbBr4_agm005739369_best_existing | 6 | 80 | 30 min |
| 04_LiBr_rs_agm003157565_ref | 2 | 80 | 5 min |
| 05_SbBr3_agm002139075_ref | 8 | 80 | 35 min |
| 06_LiBiBr4_TP01_hit | 12 | 80 | 90 min |
| 07_LiBiBr4_agm005739442_best_existing | 6 | 80 | 30 min |
| 08_BiBr3_agm002138960_ref | 8 | 80 | 35 min |

Reused, not recomputed: NaCl, SbCl3 (both polymorphs), BiCl3 at 60 Ry (night2 jobs 07-10); NaSbCl4 hit and best existing (night2 05, 06). Not needed for the NaSbCl4 verdict: NaSbCl6 (Sb(V), 16 atoms, hull 0; the balanced decomposition route is about +150 meV/atom uphill per skeptic3). LiBiCl4 completes with night2 jobs 11-15.

Overflow, inputs written to `overnight/night3/overflow/` (not auto-run; about 5.0 h more; move into `jobs/` to run): 09 NaSbBr4 hit (90), 10 NaSbBr4 agm003733885 (30), 11 NaBr agm003157587 (5), 12 Na3SbBr6 agm016963778 (60; hull 0, the competing ternary), 13 TP-LiSbCl4-03 (70), 14 LiSbCl4 agm005453208 (25), 15 SbCl3 agm002142958 at 80 Ry (20; needed because the night2 SbCl3 is at 60 Ry and LiCl is at 80 Ry). LiCl 80 Ry comes from night2 job 12.

Br pseudopotential: `Br.pbe-n-kjpaw_psl.1.0.0.UPF`, MD5 b00de77caabc69c0a6b76eb37264f235, in `overnight/pseudo/`; added by loading `make_inputs.py` source with Br in the species tuple (file not edited). No Br job has been run before, so job 02 or 04 is also the first test of that pseudopotential.

## Needs pseudopotential (nothing downloaded)

Proposed PSlibrary 1.0.0 PBE PAW file names (confirm against the PSlibrary listing before downloading):

| Element | Proposed file | Blocks |
|---|---|---|
| Al | Al.pbe-n-kjpaw_psl.1.0.0.UPF | LiAlBr4 (kept hit, ML decomposition -48) |
| Ga | Ga.pbe-dn-kjpaw_psl.1.0.0.UPF | NaGaI4, Ga chlorides/bromides |
| In | In.pbe-dn-kjpaw_psl.1.0.0.UPF | LiInBr4, LiInCl4 |
| Sc | Sc.pbe-spn-kjpaw_psl.1.0.0.UPF | LiScBr4, LiScCl4 |
| Y | Y.pbe-spn-kjpaw_psl.1.0.0.UPF | LiYBr4 |
| I | I.pbe-n-kjpaw_psl.1.0.0.UPF | NaGaI4 and all iodides |

Candidate files: `discovery/candidates/<ID>/gate3_check.json` and `relaxed_mace_d3_rattled_gate3.cif` added for each kept hit; `relaxed_mace_d3.cif` and `provenance.json` were already present.
