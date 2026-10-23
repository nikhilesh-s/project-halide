# NaBiCl4: adversarial search for anything below the KBiF4-type winner (cycle 5)

Everything here is a **PROPOSED MODEL** ranked by a machine-learned potential (MACE-MPA-0 medium, + D3(BJ) float64 for the final numbers). No DFT was run in this cycle; nothing is a stability or discovery claim.

Reference = the DFT-tested winner (KBiF4-type transplant of agm077914577, DFT-relaxed geometry from `overnight/night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.out`, re-relaxed with MACE+D3 after a 0.03 A rattle). All dE are meV/atom relative to it.

## Calibration (how far to trust the ML ranking)

- MACE+D3: winner - agm088474086 = **-5.47 meV/atom**; DFT (E31, 80 Ry): **-5.71**.
- MACE+D3: winner - agm084661913 (its Fdd2 minimum) = -2.09; DFT: -9.56. MACE+D3: winner - agm005743614 (C2/c) = -18.14; DFT: -25.84.
- So the model reproduces one gap to 0.2 meV/atom and misses the two others by about 7.5; skeptic cycle 3 measured an RMS of 3.4 meV/atom for Na-Bi-Cl. **Any ML margin below roughly 5-8 meV/atom is not a ranking**; it is a reason to run DFT.

## Reference set at the same level

| structure | atoms | dE MACE+D3 | V A^3/atom | SG (0.1) | dE no-D3 cheap | Na-Cl (<= 3.3 A) | Bi-Cl (<= 3.2 A) |
|---|---|---|---|---|---|---|---|
| RIV_KGaI4type | 24 | -2.64 | 26.12 | P2_12_12_1 | +10.68 | CN7 (2.76-3.24); CN7 (2.77-3.23); CN7 (2.76-3.25); CN7 (2.76-3.26) | CN6 (2.61-2.95) |
| WIN_KBiF4type_MLgeom | 12 | -0.01 | 28.19 | P2_12_12 | +0.15 | CN6 (2.74-2.87); CN6 (2.75-2.87) | CN6 (2.62-2.88); CN6 (2.62-2.87) |
| WIN_KBiF4type_DFTgeom | 12 | +0.00 | 28.16 | P2_12_12 | +0.00 | CN6 (2.75-2.87); CN6 (2.75-2.86) | CN6 (2.62-2.87) |
| RIV_NaGaF4type | 24 | +1.87 | 28.47 | Pbcn | -0.69 | CN6 (2.76-2.86); CN6 (2.77-2.86) | CN6 (2.61-2.89) |
| agm084661913_sk3geom | 12 | +2.09 | 28.08 | Fdd2 | +3.32 | CN6 (2.76-2.84) | CN6 (2.61-2.86) |
| agm088474086 | 12 | +5.47 | 28.00 | P2/c | +8.03 | CN6 (2.77-2.82) | CN6 (2.61-2.87) |
| agm084661913 | 12 | +7.75 | 28.86 | Fdd2 | +3.53 | CN6 (2.74-2.86); CN6 (2.74-2.85) | CN6 (2.61-2.86) |
| agm005743614 | 12 | +18.14 | 30.49 | C2/c | -3.66 | CN6 (2.80-2.83) | CN6 (2.61-2.89) |
| agm003736714 | 6 | +25.17 | 28.43 | P2/m | +17.82 | CN6 (2.78-2.89) | CN6 (2.69-2.78) |
| agm075995748 | 6 | +26.41 | 28.74 | Cmmm | +21.81 | CN6 (2.74-2.83) | CN6 (2.67-2.77) |

The KGaI4-type rival (24 atoms, already queued for DFT) is **below** the winner at this level; agm084661913 has two minima depending on the rattle (the lower one is the skeptic-cycle-3 geometry).

## Throughput caveat

`taskpolicy -b` (required: the laptop was running DFT and in use) puts the workers on throttled efficiency cores: a relaxation step cost about 10x what it did in the cycle-3 scan. The counts below are therefore far smaller than the few hundred to 1500 random starts that were asked for. The scripts are resumable: rerunning them with longer deadlines extends every stage.

## Strategy 1: wider prototype pool

- Sources: Alexandria PBE (one anonymous-formula A4BC query + the cycle-3 cache) and Materials Project; ABX4 with A in Li/Na/K/Rb/Cs/Ag/Tl/Cu, B in Al/Ga/In/Tl/Sc/Y/La-Lu/Bi/Sb/As/Fe/Cr/Au, X in F/Cl/Br/I, <= 48 atoms, hull <= 60 meV/atom: {'alexandria': 2327, 'mp': 192} entries. OQMD: unreachable: <HTTPError 502: 'Bad Gateway'>.
- After role-preserving dedup (StructureMatcher default, inside equal primitive size + space group buckets): 1006 duplicates of the 228 prototypes already tried, 431 duplicates among the rest, **359 new prototypes** (by primitive atoms: {6: 72, 12: 138, 18: 23, 24: 68, 36: 47, 48: 11}).
- Cheap pass (MACE float32, no D3, LBFGS, fmax 0.05, <= 150 steps, start volume 28 A^3/atom), <= 24 atoms and low hull first until the deadline: **110 of 359 relaxed** ({'alexandria': 106, 'mp': 4}; by atoms {6: 26, 12: 71, 18: 13}); failed/collapsed 0. **249 new prototypes were NOT tried** (all larger cells / higher hull).
- Re-found known NaBiCl4 structures: {'agm075995748': 7, 'agm084661913': 3, 'agm003736714': 2, 'WIN_KBiF4type': 1, 'agm005743614': 1}.

## Strategy 2: random structure search (AIRSS-style)

- Starts generated: 240 (Z: {1: 61, 2: 125, 4: 54}; 144 with imposed symmetry; 22 needed the distance constraints relaxed by 5 %). Z = 3 never came out of the generator; Z = 4 is under-sampled.
- Relaxed (MACE float32, no D3, LBFGS, fmax 0.05, <= 250 steps): **33** (Z: {1: 7, 2: 13, 4: 13}); usable 33, converged 25.
- Distinct minima among the usable ones (same primitive size, |dE| < 4 meV/atom, StructureMatcher): **33**; 0 of them were hit more than once.
- **Rediscovery rate**: winner re-found by 0 of 33 starts; known structures re-found at least once: 3 of 8 ({'agm075995748': 1, 'agm088474086': 2, 'agm003736714': 1}). Known set = winner, KGaI4-type and NaGaF4-type rivals, five Alexandria entries.
- Reading: a search that does not re-find the known low-energy structures, and whose minima are almost all hit once, has not sampled this landscape; its failure to find something lower is weak evidence.

## Strategy 3: local perturbations (supercell + rattle + strain)

| job | parent | atoms | level | dE vs parent (meV/atom) | converged | outcome |
|---|---|---|---|---|---|---|
| pt_WIN_211_r0.1_p3 | WIN_KBiF4type_DFTgeom | 24 | cheap no-D3 | +0.01 | True | returns to parent |
| pt_WIN_211_r0.3_m3 | WIN_KBiF4type_DFTgeom | 24 | cheap no-D3 | +0.03 | True | returns to parent |
| pt_WIN_211_r0.2_d3 | WIN_KBiF4type_DFTgeom | 24 | MACE+D3 | -0.01 | True | returns to parent |
| pt_WIN_121_r0.1_p3 | WIN_KBiF4type_DFTgeom | 24 | cheap no-D3 | +0.05 | True | returns to parent |
| pt_WIN_121_r0.3_m3 | WIN_KBiF4type_DFTgeom | 24 | cheap no-D3 | +0.04 | True | returns to parent |
| pt_WIN_121_r0.2_d3 | WIN_KBiF4type_DFTgeom | 24 | MACE+D3 | -0.01 | True | returns to parent |
| pt_WIN_112_r0.1_p3 | WIN_KBiF4type_DFTgeom | 24 | cheap no-D3 | +0.03 | True | returns to parent |
| pt_WIN_112_r0.3_m3 | WIN_KBiF4type_DFTgeom | 24 | cheap no-D3 | -0.01 | True | returns to parent |
| pt_WIN_112_r0.2_d3 | WIN_KBiF4type_DFTgeom | 24 | MACE+D3 | -0.01 | True | returns to parent |
| pt_KGaI4_111_r0.1_p3 | RIV_KGaI4type | 24 | cheap no-D3 | -0.48 | True | returns to parent |
| pt_KGaI4_111_r0.3_m3 | RIV_KGaI4type | 24 | cheap no-D3 | -0.40 | True | returns to parent |
| pt_KGaI4_211_r0.2_p3 | RIV_KGaI4type | 48 | cheap no-D3 | -0.41 | True | returns to parent |
| pt_KGaI4_121_r0.2_m3 | RIV_KGaI4type | 48 | cheap no-D3 | -0.39 | True | returns to parent |
| pt_KGaI4_112_r0.2_p3 | RIV_KGaI4type | - | cheap | - | - | not run (deadline) or failed |
| pt_NaGaF4_111_r0.1_p3 | RIV_NaGaF4type | - | cheap | - | - | not run (deadline) or failed |
| pt_NaGaF4_111_r0.3_m3 | RIV_NaGaF4type | - | cheap | - | - | not run (deadline) or failed |
| pt_NaGaF4_211_r0.2_p3 | RIV_NaGaF4type | - | cheap | - | - | not run (deadline) or failed |
| pt_NaGaF4_121_r0.2_m3 | RIV_NaGaF4type | - | cheap | - | - | not run (deadline) or failed |
| pt_NaGaF4_112_r0.2_p3 | RIV_NaGaF4type | - | cheap | - | - | not run (deadline) or failed |

## MACE+D3 stage

- Triage score = dE_cheap + k (V_cheap - V_cheap,winner), k = 5.1 meV/A^3 fitted on the reference set (rms 2.6 meV/atom). Reason: D3 favours dense cells; the no-D3 energy alone ranks the KGaI4-type rival at +10.7.
- Cheap pool 153; inside the 25 meV/atom score window (widened from 15 because so few candidates existed) 26; distinct 10; already-known structures 5; sent to MACE+D3 (float64, LBFGS + cell, fmax 0.01, 0.03 A rattle): 5; **finished before the deadline: 5**.

### Distinct structures within +15 meV/atom of the winner (MACE+D3), lowest first

| # | dE meV/atom | atoms | V | SG 0.01 / 0.1 | identity vs winner, rivals, 5 Alexandria | Na-Cl | Bi-Cl | conv | phonon 2x2x2 | provenance |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | -0.01 | 12 | 28.19 | P1 / P2_12_12 | WIN_KBiF4type_DFTgeom(default) | CN6 (2.75-2.87) | CN6 (2.62-2.87) | True | not run / not finished | S3 perturbation pt_WIN_121_r0.2_d3 (+2 more) |
| 2 | +0.65 | 12 | 25.85 | P-1 / P-1 | no match (new) | CN7 (2.80-3.01); CN7 (2.81-3.01) | CN7 (2.62-2.92) | True | min -0.00 THz at q=[0.0, 0.0, 0.0] (0 modes < -0.3) | S1 transplant of NaHoCl4 mp-aaacqvxg (mp, 12 atoms, hull 0 meV) |
| 3 | +7.00 | 24 | 28.77 | P-1 / Pbcn | no match (new) | CN6 (2.77-2.80) | CN6 (2.61-2.89) | True | not run / not finished | S2 random start rs_0004 (Z=4, start symmetry P1) |
| 4 | +8.35 | 12 | 28.73 | P-1 / P2/c | no match (new) | CN6 (2.80-2.83) | CN6 (2.62-2.90) | True | not run / not finished | S2 random start rs_0024 (Z=2, start symmetry P1) |

## Deliverables

- NBC-01: dE +0.65 meV/atom, 12 atoms, V 25.85, S1 transplant of NaHoCl4 mp-aaacqvxg (mp, 12 atoms, hull 0 meV); QE input in `overnight/later/nbc_rivals/`

New distinct structures at or below the winner: **0**; within +5 meV/atom: 1. Already known and still untested by DFT: KGaI4-type rival at -2.64 (24 atoms; `discovery/candidates/TP-NaBiCl4-01/`), NaGaF4-type at +1.87.

## Judgement

- Nothing new was found below the winner, but the search was thin (see rediscovery rate and the untried prototypes), so this is weak support, not confirmation.
- The real threats are the dense, Na-CN7 structures (V about 26 A^3/atom): the KGaI4-type rival (below the winner at MACE+D3) and any NBC entry listed above. Their ranking leans on D3(BJ); the QE runs use zero-damping D3, which was not tested here and may shift dense cells by several meV/atom either way.
- With margins of 0-3 meV/atom against a model error of 3-8 meV/atom, 'lowest-energy NaBiCl4 structure found' is undecided until those DFT jobs run. Defensible wording today: 'lowest of the structures tested by DFT'.

## What was NOT covered

- Z > 4 (more than 24 atoms) in the random search; Z = 3 not sampled at all; Z = 4 thinly. Random starts: tens, not hundreds.
- 249 of the 359 new prototypes (mostly 18-48 atoms) never got a cheap relaxation; OQMD was unreachable (HTTP 502); ICSD, GNoME, AFLOW not queried. Prototype chemistries outside the A(+1) B(+3) X4 lists (e.g. A2+ B2+ X4, oxides/chalcogenides ABX4) were not used as templates.
- Non-stoichiometric, disordered, partially occupied or defect models; other Na-Bi-Cl compositions (the hull question is separate: Na3BiCl6 + BiCl3 is covered by DFT in E28/E31).
- Bi lone-pair physics beyond what MACE-MPA-0 learned from PBE; no spin-orbit coupling; no finite temperature / vibrational free energy; no pressure.
- Triage was done without D3 (with a volume correction): a structure whose basin only exists with dispersion can be missed. Cheap relaxations were capped (150 / 250 steps), unconverged ones were kept in the ranking only by their last energy.
- Phonons: only the structures listed with a phonon entry; 24-atom contenders (192-atom supercells) were not attempted.
- Everything is in-sample-biased ML: the Alexandria entries are in MACE's training set, the new models are not.

## Compute

- Core time in relaxations: cheap 180 core-min, MACE+D3 31 core-min; 4 single-thread workers under `taskpolicy -b`.
- Scripts: `discovery/scripts/nbc_search/` (common, s1_fetch, s1_run, s2_airss, s3_perturb, d3_final, phonon, report, chain.sh). Caches: `discovery/cache/nbc_search/`.
