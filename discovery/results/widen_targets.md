# Widening the search: target compounds beyond Li/Na tetrahalides (2026-10-13)

Everything computed here is a **PROPOSED MODEL** ranked by a machine-learned potential (MACE-MPA-0 + D3). No DFT was run. Nothing below is a stability or discovery claim. "Not found" means not found in the listed queries and is NOT clearance: ICSD, Pauling File/MPDS, PDF-4+, Gmelin and Zh. Neorg. Khim. were not accessible.

Scripts: `discovery/scripts/widen/db_sweep.py` (Alexandria + Materials Project + MC3D for all 180 A-M-X systems), `db_table.py` (-> `discovery/results/widen/db_table.csv`, 1 row per composition), `widen_scan.py` (Part B/C). Raw relaxations: `discovery/results/widen/{cheap,d3}/`; machine summary `widen/scan_summary.json`; per-target tables `widen/partB.md` (first reduced run kept as `partB_run1_short12.md`).

## Part A: target table

Database columns: Alexandria PBE = number of listed entries / lowest hull distance (meV/atom) / atoms in that entry. MP = Materials Project summary for the chemical system. MC3D = Materials Cloud MC3D (structures imported from experimental databases; coverage is incomplete, e.g. K2BiCl5 (2021) is absent). MP's `theoretical` flag is unreliable (E32).

### Class 1: KNOWN (phase-diagram compound) + STRUCTURE UNREPORTED in searched sources

| # | compound | evidence it exists | structure search | Alexandria | ion-conductor relevance | note |
|---|---|---|---|---|---|---|
| 1 | KBiCl4 | Tabiev 1982, Arm. Khim. Zh. 35(7) (https://arar.sci.am//Content/181139/file_0.pdf, E24): treated as established, NaBiCl4-KBiCl4 section studied by DTA/XRD. Atomics International report NAA-SR-5093 (https://www.osti.gov/biblio/1963097): KCl-BiCl3 has K3BiCl6 (congruent) and K2BiCl5 (incongruent) and "another compound such as KBiCl4 may exist" (search-tool reading of an OCR scan; PDF not text-readable here). CONTESTED: Zhou et al. 2021 (doi:10.1016/j.jssc.2021.122621, snippet) say several studies agree on only K3BiCl6 and K2BiCl5 | none found: MP K-Bi-Cl has only a theoretical K3BiCl6; MC3D empty; web queries Q1, Q5, Q11 | 4 / 1.5 / 12 (C2/c) | low for K+; direct structural analogue of NaBiCl4 and of Tabiev's NaBiCl4.KBiCl4 | screened |
| 2 | RbBiCl4 | Tabiev 1982 (as above); Addison & Halstead 1966 (doi:10.1039/J19660001236, abstract): solid complex formation with RbCl, formula not seen | none found: MP and MC3D empty for RbBiCl4; Q5, Q13 | 4 / 0.0 / 24 (P2_1/c) | low | screened |
| 3 | K3BiCl6 | Dzeranova & Trunin 2007 (doi:10.1134/S002016850709018X): 3KCl.BiCl3 congruently melting; NAA-SR-5093; Zhou 2021 | no anhydrous structure found (Q11): MP entry is a theoretical Fm-3m decoration; MC3D empty. Known relatives: Cs3BiCl6 and CsK2BiCl6 (Benachenhou et al. 1986, J. Solid State Chem., https://www.sciencedirect.com/science/article/abs/pii/002245968690085X), K7Bi3Cl16 from solution (Kelly et al. 2016, J. Alloys Compd., https://www.sciencedirect.com/science/article/abs/pii/S0925838816303243) | 7 / 0.0 / 20 (R-3) | low | screened |
| 4 | LiBi4Cl13 (LiCl.4BiCl3) | Dzeranova & Trunin 2007, congruently melting. CONTESTED: Addison & Halstead 1966 report no compound with LiCl | absent from Alexandria, MP, MC3D | none listed | low (Li-poor) | 18 atoms per formula unit and no A M4 X13 prototype pool: NOT screened |
| 5 | K2BiBr5 | proposed in KBr-BiBr3 "with no structural information" (Zhou 2021, as summarised by the search tool; the primary phase-diagram paper was not identified) | none found (Q3, Q6) | 2 / 16.7 / 32 (Pnma, i.e. the K2PrCl5 type Zhou found for K2BiCl5) | low | needs 32 atoms: NOT screened. Likely isotypic with K2BiCl5; cheap DFT check of the listed entry |
| 6 | K4BiBr7 | same source as row 5 | absent from Alexandria, MP, MC3D | none listed | low | no pool: NOT screened. Could really be K7Bi3Br16-type (cf. K7Bi3Cl16, Rb7Sb3Br16); that is a guess |
| 7 | Rb3BiCl6 (anhydrous) | weak: Kelly 2016 characterise the hemihydrate Rb3BiCl6.0.5H2O; Addison & Halstead 1966 report RbCl-BiCl3 solid complexes (formula not seen) | anhydrous structure not found; MP theoretical only | 5 / 0.0 / 20 (R-3) | low | screened; class 1 only if the anhydrous phase is confirmed |
| 8 | LiBiCl4 (LiCl.BiCl3) | Dzeranova 2007 ("thermally unstable"), contested (E23, scout 3) | none found | 2 / 40.1 / 6 | Li conductor candidate | already in the pipeline (TP-LiBiCl4-01); not rescreened |

NaBiCl4 (the confirmed case, E28) belongs to this class and is not repeated.

### Class 2: NOT REPORTED as a compound in searched sources

| # | compound | what was found | Alexandria | note |
|---|---|---|---|---|
| 9 | KSbCl4, RbSbCl4, CsSbCl4 | no SbCl3-MCl phase diagram located (Q2, Q4, Q9, Q14); K2SbCl5 is experimentally known (MP, ICSD-linked, P2_1/c, 32 atoms; MC3D mpds S1411936) | 3 / 12.3 / 12; 2 / 11.1 / 24; 4 / 8.5 / 24 | KSbCl4 screened (hits, see Part B); RbSbCl4 not finished |
| 10 | KBiBr4, RbBiBr4 | no report found (Q3, Q6) | 3 / 4.2 / 12; 3 / 9.1 / 12 | KBiBr4 screened (hit) |
| 11 | CsBiCl4 | Dzeranova 2007 lists only 3CsCl.BiCl3 for Cs; Cs3BiCl6 and Cs3Bi2Cl9 are known structures | 3 / 9.7 / 12 | screened: ML puts it ABOVE the Cs3BiCl6/Cs2BiCl5 + BiCl3 line |
| 12 | CsBiBr4 | EVIDENCE AGAINST: "could not be verified experimentally" in the RbBr/CsBr-BiBr3 reinvestigation (Z. Anorg. Allg. Chem. 2016, doi:10.1002/zaac.201500820, snippet) | 4 / 20.8 / 12 | skip |
| 13 | KInCl4, RbInCl4, CsInCl4 | nothing found in one query (Q12); InCl3-MCl diagrams (Fedorov school; CALPHAD InCl3-KCl in Mater. Chem. Phys. 2015) not read. All three sit ON the Alexandria hull with no MP/MC3D entry | 5 / 0.0 / 12; 7 / 0.0 / 24; 8 / 0.0 / 24 | status UNKNOWN, not "unreported"; best next scouting target |
| 14 | NaInCl4, LiInCl4, LiScCl4, NaYCl4 | Na3InCl6 (P-31c), Li3InCl6, Li3ScCl6, Na3YCl6 are the known phases; NaYCl4 only as a theoretical MP entry | 3 / 12.6 / 6; 3 / 13.1 / 12; 3 / 23.4 / 6; 7 / 5.1 / 12 | LiInCl4 already queued; nothing new for Li/Na here |
| 15 | K3Bi2Cl9, K3Sb2Cl9, Rb3Sb2Cl9 | K3Bi2Cl9: existence "still uncertain" (Zhou 2021 snippet). Sb ones: literature not checked | 1 / 0.0 / 56 each (2 entries for the Sb ones) | 56-atom cells; needed as DFT competitors |

### Class 3: STRUCTURE KNOWN (rediscovery only)

K2BiCl5 (Pnma, K2PrCl5 type) and K3Bi2Br9 (P2_1/c, Tl3Bi2Br9 type): Zhou, Qiu, Mishra, Mar 2021. RbBi2Br7, CsBi2Br7: ZAAC 2016 (MP: experimental, ICSD-linked). Rb3BiBr6, Cs3BiBr6 polymorphs: Chang et al. 2021 (doi:10.1002/zaac.202000402). Cs3BiCl6, Cs3Bi2Cl9, Cs3Sb2Cl9, Rb3Sb2Br9, A3M2I9 (M = Sb, Bi; A = Rb, Cs). K2SbCl5. NaScCl4 (Pbcn), Na3ScCl6, Na3ScBr6, Na3YBr6, Na3InCl6, NaInBr4, NaInI4, LiInI4, and the Al/Ga tetrahalides flagged EXP in `db_table.csv`.

### Queries (all "not found" statements refer to these)

Q1 "KBiCl4 OR RbBiCl4 OR CsBiCl4 crystal structure ... BiCl3-KCl phase diagram"; Q2 "SbCl3-KCl SbCl3-CsCl SbCl3-RbCl phase diagram compounds KSbCl4 CsSbCl4 CsSb2Cl7 K2SbCl5"; Q3 "BiBr3-KBr BiBr3-NaBr BiI3-KI BiI3-NaI phase diagram ... KBiBr4 K3BiBr6 NaBiI4 Na3BiI6"; Q4 "CsSbCl4 OR KSbCl4 OR RbSbCl4 OR CsSb2Cl7 OR KSb2Cl7 ..."; Q5 "CsBiCl4 OR RbBiCl4 OR KBiCl4 OR CsBi2Cl7 OR KBi2Cl7"; Q6 "KBiBr4 OR RbBiBr4 OR CsBiBr4 OR K2BiBr5 OR KBiI4 OR NaBiI4 OR K4BiBr7"; Q7 "InCl3-NaCl InBr3-NaBr phase diagram compounds NaInCl4 Na3InCl6 NaInBr4 Na2InCl5 Fedorov"; Q8 "Lost horses on the frontier K2BiCl5 K3Bi2Br9 ..." (twice); Q9 "SbCl3 CsCl RbCl KCl system thermal analysis compounds Cs3Sb2Cl9 CsSbCl4 ..."; Q10 "InCl3-KCl InCl3-NaCl InCl3-LiCl phase diagram compounds ..."; Q11 "K3BiCl6 anhydrous crystal structure space group"; Q12 "KInCl4 OR RbInCl4 OR CsInCl4 OR NaInCl4 crystal structure"; Q13 "Addison Halstead 1966 bismuth trichloride alkali chlorides ... Rb3BiCl6 RbBiCl4 Rb2BiCl5" and "RbCl-BiCl3 CsCl-BiCl3 phase diagram compounds ..."; Q14 "SbBr3 SbI3 SbCl3 with KBr KI KCl CsCl ... Kendall Crittenden Miller 1923".
Fetch failures (403 or unreadable PDF): ResearchGate and ScienceDirect pages for Kelly 2016 and Zhou 2021, Wiley (ZAAC 2016), OSTI 4513889 (ORNL halide phase-diagram compilation, image-only PDF), OSTI 1963097, FactSage FTsalt list. Everything attributed to those papers above is snippet-level.
Not covered at all: Bi/Sb iodide phase diagrams (only hydrates surfaced), La, Zn, Mg, Sn, Pb systems beyond the database sweep, Kaloev/Korshunov originals.

## Part B: ML screen

Method = `transplant_scan.py` logic with these differences, all forced by time: (i) prototype pool = Alexandria A'M'X'4 or A'3M'X'6 entries with A' in Li..Cs, Ag, Tl; M' in Bi, Sb, Al, Ga, In, Sc, Y; X' in F..I, <= 24 atoms (NOT every Alexandria entry of the same anonymous formula); pool hull <= 30 meV for A M X4 (cycle-3 grouping reused, 228 prototypes), <= 40 meV for A3 M X6 (43 prototypes); (ii) only 40 shortlisted prototypes per family (those already realised by a Bi/Sb chloride, then lowest cycle-3 cheap dE, then lowest hull); (iii) only the 2 lowest cheap-pass survivors per target got MACE+D3, plus the 3 lowest listed entries after a 0.03 A rattle; (iv) competitors = lowest-hull AX and MX3 plus listed ternaries on the AX-MX3 line within 5 meV of the Alexandria hull and <= 40 atoms; the 56-atom A3M2X9 entries were NOT relaxed, so "vs all competitors" is incomplete wherever one is listed. Energies in meV/atom; negative decomposition = below the competitors. Compute used: about 1 h 50 min wall (2 workers, then 4 after `nbc_search` ended).

| target (class) | best listed, MACE+D3 (SG, atoms, V A^3/atom) | Alexandria's own lowest entry vs that | best transplant (prototype) | dE vs best listed | vs binaries | vs relaxed competitors | identity | verdict |
|---|---|---|---|---|---|---|---|---|
| KBiCl4 (1) | agm082986275 (Pnna, 24, 31.09) | agm005744494 +21.2 | TlSbCl4 type, P2_1/c, 24, V 31.20 | +0.0 | -37.1 | -9.9 (K3Bi2Cl9 missing) | no match at default/loose tolerance, but degenerate in energy | no hit; listed structure stands |
| RbBiCl4 (1) | agm003739447 (P-1, 12, 32.48) | agm089761403 +7.2 | TlSbCl4 type, P2_1/c, 24, V 33.54 | +8.4 | -49.7 | -13.1 | no match | no hit |
| K3BiCl6 (1) | agm003569199 (P2_1/c, 20, 31.59) | agm058388393 +15.4 | Na3InCl6 type, P-31c, 20, V 34.64 | +18.6 | -30.4 | -8.1 (K3Bi2Cl9 missing) | no match | no hit |
| Rb3BiCl6 (1?) | agm003572728 (P2_1/c, 20, 34.94) | agm048996295 +19.8 | Rb3InCl6 type, R-3c, 20, V 35.08 | +0.6 | -65.1 | -30.2 | no match | no hit |
| CsBiCl4 (2) | agm003739771 (C2/c, 12, 36.08) | same | CsSbCl4 type, Cmc2_1, 24, V 35.40 | -6.0 | -68.9 | +2.5 (above Cs3BiCl6/Cs2BiCl5 + BiCl3) | no match | fails competitors |
| KSbCl4 (2) | agm003737289 (C2/c, 12, 30.30) | same | LiGaI4 type, P2_1/c, 24, V 30.74 = **W-KSbCl4-01**; RbBiCl4 type, P2_1/c, 24 = **W-KSbCl4-02** (-7.5 / -41.4 / -27.9) | -12.2 | -46.0 | -32.6 (K3Sb2Cl9 and K2SbCl5 missing) | no match | ML hit, PROPOSED MODEL |
| KBiBr4 (2) | agm003742469 (C2/c, 12, 35.89) | same | RbBiCl4 type, P2_1/c, 24, V 36.19 = **W-KBiBr4-01** | -10.1 | -14.9 | -13.6 (K3Bi2Br9 missing) | no match | ML hit, PROPOSED MODEL |
| RbSbCl4 (2) | not finished (cheap-pass deadline) | | | | | | | rerun |

Readings:
- For all four class-1 targets no transplant beats the listed structures. The useful result there is a re-ranking: MACE+D3 puts a different listed entry 7-21 meV/atom below the one Alexandria ranks lowest. All are in the model's training set, so this needs DFT.
- All class-1 chlorides are strongly negative against the binaries at ML level (-37 to -66), consistent with them being real compounds.
- The three hits are in class 2 (composition not reported), and their competitor sets are incomplete: K2SbCl5 is an experimentally known phase that lies 12 meV above the Alexandria hull and was therefore outside the 5 meV rule, and the 56-atom K3Sb2Cl9 / K3Bi2Br9 (the latter experimentally known) were not relaxed. Until those are included, "below all competitors" is not established even at ML level.
- Caveats of cycle 3 apply unchanged (in-sample potential, BJ-damped D3, no phonons, cheap-pass triage without D3).

## Part C: hand-off

Candidate folders (CIF + provenance.json, label PROPOSED MODEL): `discovery/candidates/W-KSbCl4-01`, `W-KSbCl4-02`, `W-KBiBr4-01`.

**No QE inputs were written** (`overnight/later/widen/` is empty): every hit contains K, and no K, Rb or Cs pseudopotential is registered in `overnight/pseudo_md5.txt`. Expected PSlibrary 1.0.0 names, to be verified and md5-registered first: `K.pbe-spn-kjpaw_psl.1.0.0.UPF`, `Rb.pbe-spn-kjpaw_psl.1.0.0.UPF`, `Cs.pbe-spn-kjpaw_psl.1.0.0.UPF`. `overnight/make_inputs.py` also needs K/Rb/Cs/Br added to its MASS, UPF and element-order lists.

DFT jobs needed (vc-relax, PBE-D3, same settings as night 2; chlorides 60 Ry, bromide 80 Ry). Atom counts are primitive cells.

| group | job | atoms | pseudopotentials (flag = not registered) |
|---|---|---|---|
| KSbCl4 | W-KSbCl4-01 | 24 | K (MISSING), Sb, Cl |
| | W-KSbCl4-02 | 24 | same |
| | best listed agm003737289 | 12 | same |
| | KCl agm003269520 | 2 | K (MISSING), Cl |
| | SbCl3 agm002205242 (already run in night 2, job 07) | 16 | Sb, Cl |
| | K3SbCl6 agm050704424 | 10 or 20 (check entry) | K (MISSING), Sb, Cl |
| | K2SbCl5, experimental structure (MP mp-aaaaboht / Alexandria K2SbCl5 P2_1/c) | 32 | same |
| | K3Sb2Cl9 agm012094232 | 56 | same |
| KBiBr4 (80 Ry) | W-KBiBr4-01 | 24 | K (MISSING), Bi, Br |
| | best listed agm003742469 | 12 | same |
| | KBr agm003157738 | 2 | K (MISSING), Br |
| | BiBr3 agm002138960 | see entry | Bi, Br |
| | K3BiBr6 agm058292496 | 10 or 20 (check entry) | K (MISSING), Bi, Br |
| | K3Bi2Br9 agm011997501 (experimentally known type) | 56 | same |
| KBiCl4 (class 1, no transplant hit) | agm082986275 (Pnna) | 24 | K (MISSING), Bi, Cl |
| | TlSbCl4-type transplant `widen/d3/tp_KBiCl4_agm005759109.json` | 24 | same |
| | agm005744494 (Alexandria lowest) | 12 | same |
| | K3BiCl6 agm003569199 | 20 | same |
| | BiCl3 agm002142654 (done, night 2), KCl | 16, 2 | |
| RbBiCl4 (class 1) | agm003739447 (P-1) | 12 | Rb (MISSING), Bi, Cl |
| | agm089761403 (P2_1/c) | 24 | same |
| | RbCl agm003157439, Rb3BiCl6 agm003572728 | 2, 20 | same |

## Which targets are worth DFT time

1. **KBiCl4** (class 1). It is the only widened target that is both reported in the literature and a direct analogue of the confirmed NaBiCl4 case. The ML screen produced no new model, but it says the lowest-energy listed structure is not the one Alexandria ranks lowest (Pnna 24-atom vs C2/c 12-atom, 21 meV/atom apart), and a TlSbCl4-type cell is degenerate with it. Four small jobs settle which listed structure is the defensible "first structural model". Caveat: its existence is contested in the literature.
2. **RbBiCl4** (class 1). Same argument, two jobs (12 and 24 atoms) plus references; costs one more pseudopotential (Rb).
3. **KSbCl4 (W-KSbCl4-01)** (class 2). Largest ML margin of this cycle (-12.2 vs best listed), but the composition is unreported and the experimentally known K2SbCl5 plus K3Sb2Cl9 must be in the DFT hull before anything is said; without them the result is uninterpretable.

Lower priority: KBiBr4 (hit, but only -14.9 vs binaries and K3Bi2Br9 untested), K3BiCl6/Rb3BiCl6 (listed structures stand; 20-atom re-ranking jobs only), CsBiCl4 (fails at ML level), LiBi4Cl13 and K4BiBr7 (need a structure search, not a transplant scan).
Not worth DFT now: anything in class 3.
Cheapest unlock: registering the K pseudopotential enables rows 1 and 3 and KBiBr4 at once.
None of this adds a Li or Na conductor candidate: the widened top-class targets are K/Rb compounds, valuable as structure-model results, not as electrolytes.
Before any claim: ICSD / Pauling File check for KBiCl4, RbBiCl4, K3BiCl6 (via the university library), and the Tabiev 1982 main text.
