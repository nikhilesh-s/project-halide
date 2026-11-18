# New candidate compositions for the last compute window (2026-10-15)

Every TP-* structure here is a **PROPOSED MODEL** from the cycle-3 prototype-transplant scan, ranked by MACE-MPA-0 + D3(BJ). No DFT has been run on any of them. Nothing here is a stability or discovery claim.

## Bottom line

The 256 scan hits contain only 20 compositions. Twelve of them already have a candidate in the DFT queue, or are dead (LiBiCl4, LiSbCl4) or known structures (LiAlBr4, LiInBr4). That leaves **8 new compositions**, and **only one has a realistic G3 chance: NaBiBr4** (ML decomposition -0.9 meV/atom, borderline). The other seven are above the AX + MX3 line at ML level by +1.4 to +18.6 meV/atom. One of 256 is not the "3-5" the brief hoped for. A new scan (Li iodides and Na-M-I were dropped for budget in cycle 3) is the only way to get more, and there is no time for that before the 2026-10-21 cutoff.

What was queued (front of `overnight/m5_oct16/jobs`, all vc-relax, PBE-D3, 80/640 Ry, k-spacing 0.22, starting from the rattled MACE+D3 cells):

| Job | What | Atoms | Est. M5 min |
|---|---|---|---|
| 00a_NaBiBr4_TP01_CsAlCl4type_hit | **TP-NaBiBr4-01**, parent agm083138585 (CsAlCl4-type, P-1) | 12 | 90 |
| 00b_NaBiBr4_agm069058740_best_existing | best listed NaBiBr4 (G2) | 6 | 30 |
| 00c_NaBiBr4_agm005739654_listed2 | second listed NaBiBr4, +3.1 in ML (inside the ML error) | 6 | 30 |
| 00d_Na3BiBr6_agm053446239_ternary | Na3BiBr6, **on the Alexandria PBE hull** (G3) | 10 | 50 |
| 00e_NaSbI4_TP01_CsAlCl4type_hit | TP-NaSbI4-01, parent agm083138585 (CsAlCl4-type), long shot | 12 | 100 |
| 00f_NaSbI4_agm005743434_best_existing | best listed NaSbI4 (G2) | 6 | 35 |
| 00g_SbI3_agm003264389_ref | SbI3, lowest in MACE+D3 (G3) | 8 | 40 |

About 6.3 M5 hours (Na-Bi-Br set 3.3 h, Na-Sb-I set 2.9 h). Already computed, not repeated: NaBr agm003157587 80 Ry (`air/15`), R-3 BiBr3 agm002138960 80 Ry (`air/14`); C2/m BiBr3 (experimental) runs tonight (`m5_oct15/10i`); NaI agm003157588 80 Ry is `later/08` on the Air.
Not auto-run (`overnight/m5_oct16/overflow`, move into `jobs/` to run, about 2.9 h): 00h TP-NaBiI4-01 (TlAlCl4-type, 12-atom primitive), 00i NaBiI4 agm003740643, 00j BiI3 agm003263390.
`overnight/make_nights6plus.py` now lists the Bi-Br-Na and I-Na-Sb systems, so the next morning run of it writes their 80 Ry k0.17 G4 SCFs once the relaxations finish.

The queued jobs push the rest of m5_oct16 (cheap SCFs) and m5_oct17 (40-atom Li3SbBr6) back by about 6 h. If the Li3SbBr6 result matters more, move 00e-00g into `overflow/`.

## G1: rattle check (MACE-MPA-0 + D3, rattle 0.03 Å seed 3, FIRE + FrechetCellFilter fmax 0.01)

Scripts `discovery/scripts/oct15_new/{run,analyze,make_queue}.py` (copies of gate4's, with ternaries up to 30 meV/atom above the Alexandria hull so Na3BiBr6 is included). Raw: `discovery/results/oct15_new/`. 26 relaxations, all converged, under 1 min wall at `nice -n 10`.

| Hit | Atoms | dE vs best listed, scan -> after rattle | Decomposition vs AX + MX3, scan -> after | Full ML hull (incl. ternaries) | Best listed after rattle (others, meV above) | Gate-3 rule (dE <= -5 and decomp <= -3) |
|---|---|---|---|---|---|---|
| TP-NaBiBr4-01 | 12 | -17.3 -> -17.3 | -0.9 -> -0.9 | -0.9 (NaBr + BiBr3; Na3BiBr6 agm053446239 and agm064728767 are above that line in ML) | agm069058740 (+3.1 agm005739654, +12.2 agm003733888) | structure stable on rattle; decomposition fails the -3 cut |
| TP-NaSbI4-01 | 12 | -20.5 -> -20.4 | +6.2 -> +7.0 | +7.0 (NaI + SbI3) | agm005743434 (+5.1) | fails decomposition |
| TP-NaBiI4-01 | 12 | -9.6 -> -9.6 | +5.7 -> +8.0 | +8.0 (NaI + BiI3) | agm003740643 (+6.1) | fails decomposition |

All three survive the rattle as structures (no slide to another minimum; dE unchanged). None passes the gate-3 decomposition cut. They are queued anyway because this is the best the scan has, and ML decomposition energies in the Na systems have been off by about ±6 meV/atom against DFT in both directions:

| System | ML decomposition (after rattle) | DFT vs AX + MX3 | DFT minus ML |
|---|---|---|---|
| NaBiCl4 (KBiF4-type) | -13.8 | -19.8 | -6.0 |
| NaSbBr4 | -5.2 | -6.2 | -1.0 |
| NaSbCl4 | -18.1 | -12.0 | +6.1 |
| Li systems (LiBiCl4, LiSbBr4, LiBiBr4, LiSbCl4-03), for contrast | -6.9 / -19.8 / -11.1 / -5.9 | +7.7 / -8.0 / -3.7 / +15.4 | +7 to +21 (always worse) |

## Ranking of everything considered (best hit per composition, top 15 of 20)

Ranking keys in order: new composition (rule 3), not reported experimentally (rule 2), ML decomposition vs AX + MX3 (the G3 proxy; rule 4), references already done at 80 Ry (rule 1), cell size (rule 5). Scan numbers (MACE+D3, meV/atom). "Decomp" = E(AMX4) - [2 AX + 4 MX3]/6 per atom, negative = stable vs binaries.

| Rank | Composition | Best hit (parent prototype) | Atoms | dE no-D3 | dE D3 | Decomp | Novelty (scout) | Refs at 80 Ry | Decision |
|---|---|---|---|---|---|---|---|---|---|
| 1 | NaBiBr4 | TP-NaBiBr4-01, agm083138585 (CsAlCl4) | 12 | -15.0 | -17.3 | -0.9 | compound not reported (scout 4, sec. 5) | NaBr, R-3 BiBr3 done; C2/m BiBr3 tonight; Na3BiBr6 + 2 listed queued | **queued 00a-00d** |
| 2 | NaSbI4 | TP-NaSbI4-01, agm083138585 (CsAlCl4) | 12 | -16.1 | -20.5 | +6.2 (+7.0 rattled) | **not scouted** | NaI on the Air (later/08); SbI3 + listed queued | **queued 00e-00g** (long shot) |
| 3 | NaBiI4 | TP-NaBiI4-01, agm078380533 (TlAlCl4) | 24 (12 prim.) | -8.7 | -9.6 | +5.7 (+8.0 rattled) | **not scouted** | NaI on the Air; BiI3 + listed written | overflow only |
| 4 | LiScCl4 | TP-LiScCl4-01, agm089761403 (RbBiCl4) | 24 | -4.7 | -13.7 | -3.5 | not scouted; Li3ScCl6 is a known conductor and must be beaten | none (Sc chloride, Li3ScCl6 not computed) | rejected: 12-atom TP-02 is -0.8 after rattle (gate 4), TP-01 is 24 atoms, competitor set large |
| 5 | NaInCl4 | TP-NaInCl4-01, agm082661748 (NaYBr4) | 24 | -6.1 | -6.2 | +1.4 | not reported, weak evidence against (scout 4, sec. 11) | InCl3 on the Air; Na3InCl6 (known, Alexandria hull 0) not computed | rejected |
| 6 | NaScBr4 | TP-NaScBr4-01, agm082661748 (NaYBr4) | 24 | -12.9 | -8.7 | +12.2 | not scouted | none | rejected |
| 7 | NaInBr4 | TP-NaInBr4-01, agm005739607 (NaSbBr4) | 6 | +8.8 | -13.9 | +18.6 | not scouted | none | rejected |
| 8 | LiInBr4 | TP-LiInBr4-01, agm082912513 (RbGaF4) | 12 | +11.7 | -7.5 | -11.4 | **structure known** (spinel) | - | rediscovery only |
| 9 | LiAlBr4 | TP-LiAlBr4-01, agm064525236 (TlInCl4) | 12 | +3.3 | -16.3 | -47.9 | **structure known** (P2_1/c) | - | rediscovery only |
| 10 | LiSbCl4 | TP-LiSbCl4-01, agm003231612 (NaSbF4) | 24 | -54.4 | -73.5 | -15.0 | not reported | done | same composition as TP-LiSbCl4-03, which fails G3 by +15.4; with that ML-to-DFT shift (+21) TP-01 lands near +7 |
| 11 | LiBiCl4 | TP-LiBiCl4-01, agm084160006 (RbAlCl4) | 24 | -26.4 | -33.6 | -13.6 | contested | done | already failed G3 (+7.7) |
| 12 | LiSbBr4 | TP-LiSbBr4-01, agm078380533 (TlAlCl4) | 24 | -23.6 | -23.7 | -19.7 | not reported | | existing candidate (DFT model is the 12-atom one) |
| 13 | LiBiBr4 | TP-LiBiBr4-01, agm082912513 (RbGaF4) | 12 | -14.9 | -23.4 | -11.1 | not reported | | existing candidate |
| 14 | NaSbBr4 | TP-NaSbBr4-01, agm078380533 (TlAlCl4) | 24 | -22.1 | -22.0 | -5.1 | not reported | | existing candidate |
| 15 | NaGaI4 | TP-NaGaI4-01, agm003232905 (LiGaI4) | 24 | +6.1 | -21.3 | -8.1 | not scouted | | existing candidate (TP-02 queued on the Air) |

Not listed (existing candidates): NaBiCl4, NaSbCl4, LiInCl4, LiScBr4, LiYBr4.

## Why each pick

- **TP-NaBiBr4-01.** The only new composition with an ML decomposition below zero. Large ML margin over every listed NaBiBr4 (-17.3 with D3, -15.0 without, so not a D3 artefact). Small cell. Same CsAlCl4-type parent as TP-NaSbCl4-01, which DFT confirmed against its listed structures. Two of its references (NaBr, R-3 BiBr3) are already done at 80 Ry. Risks: the margin against NaBr + BiBr3 is -0.9 meV/atom, below ML resolution; Na3BiBr6 agm053446239 sits on the Alexandria PBE hull and could cut under it in DFT; experimental BiBr3 (C2/m, tonight) may lower the hull. History: NaBiBr4 on the NaBiCl4 prototype was demoted in cycle 2 (E19, H3). This is a different, lower structure (scan: no match to any Alexandria entry).
- **TP-NaSbI4-01.** Biggest ML margin over listed structures among the remaining new compositions (-20.4), 12 atoms, same CsAlCl4 parent. However it is +7.0 above NaI + SbI3 at ML level, so it needs DFT to come out about 7 meV/atom more favourable than ML. That has happened once in the Na systems (NaBiCl4, -6.0). Novelty is unscouted: Na-Sb-I was outside scout 4, so a pass would still need a literature check.
- **TP-NaBiI4-01 (overflow).** Same story as NaSbI4 with a smaller margin over listed structures (-9.6) and +8.0 vs binaries. Written but not auto-run.

## Honest odds of passing G1-G4 by 2026-10-21

The M5 should reach 00a-00g on the night of Oct 16 if night2's rival and m5_oct15 finish tonight. The G4 SCFs follow a night later.

| Candidate | G1 | G2 | G3 | G4 (convergence + scout) | All four by Oct 21 |
|---|---|---|---|---|---|
| TP-NaBiBr4-01 | pass (ML) | about 85% | about 30% | about 85% if G3 passes (needs a manual ICSD / Pauling File check, same as NaBiCl4) | **about 20-25%** |
| TP-NaSbI4-01 | pass (ML) | about 80% | about 10% | about 70% if G3 passes (unscouted) | **about 5%** |
| TP-NaBiI4-01 | pass (ML) | about 70% | about 10% | not queued | about 0% unless moved into `jobs/` |

Expected new passes from this batch: about 0.3. Plan the count of five without these.

Candidate files: `discovery/candidates/TP-{NaBiBr4,NaSbI4,NaBiI4}-01/relaxed_mace_d3_rattled_oct15.cif` (the DFT starting cell) and `oct15_check.json`; the scan's `provenance.json` is already in each folder, with parent prototype agm083138585 (CsAlCl4) for NaBiBr4 and NaSbI4 and agm078380533 (TlAlCl4) for NaBiI4. Folder numbers match the cycle-3 scan (`NAMING.md`); job names use TP01 = that folder.
