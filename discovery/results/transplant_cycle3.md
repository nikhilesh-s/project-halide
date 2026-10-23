# Cycle 3: prototype-transplant scan, A M X4 (A = Li, Na; M = Al, Ga, In, Sc, Y, Sb, Bi; X = Cl, Br, I)

Every structure here is a **PROPOSED MODEL** ranked by a machine-learned potential (MACE-MPA-0, with D3 for the final ranking). Nothing below is a stability or discovery claim; no DFT was run.

## Counts

- Alexandria entries downloaded (<= 24 atoms): 999 across 196 A'M'X'4 chemistries.
- Prototype pool (hull <= 30 meV/atom): 723 entries -> **228 distinct prototypes** (role-preserving StructureMatcher at default tolerance: A', M', X' are mapped to fixed roles before matching, which is stricter than fit_anonymous).
- Chemistry x prototype pairs: 9576; already listed for that chemistry: 145; **transplants attempted (cheap MACE pass): 9431**.
- Kept for MACE+D3 (within 15 meV/atom of best existing in the cheap pass, not matching an existing entry, de-duplicated): 1464.
- Row status counts: HIT_proposed_model: 256, already_listed: 145, cheap_above_window: 6568, dropped_budget: 907, dropped_budget_D3: 175, duplicate_of_other_transplant: 323, hit_duplicate_of_other_transplant: 96, rediscovery_D3: 35, rediscovery_cheap: 169, survives_identity: 902.
- Dropped for budget (chemistries not finished before the stage deadline): LiAlI4, LiBiI4, LiGaI4, LiInI4, LiSbI4, LiScI4, LiYI4, NaGaI4, NaInI4, NaScI4, NaYI4.
- Prototype scope this run: 228 of 228 prototypes (<= 24 atoms); 0 chemistry x prototype pairs skipped for that reason.

## Hits: survive identity check and are >= 3 meV/atom below the best existing entry at MACE+D3

dE = E(transplant) - E(lowest Alexandria entry of the same composition), both MACE+D3 relaxed; for hit chemistries every existing entry (<= 24 atoms) was relaxed with D3, not only the three lowest. E_decomp = E(AMX4) - [2 E(AX) + 4 E(MX3)]/6 per atom, positive = unstable against the binaries.

| candidate | composition | prototype source | dE no-D3 cheap (meV/atom) | dE MACE+D3 (meV/atom) | best existing (D3) | E_decomp D3 (meV/atom) | binaries AX + MX3 | SG | V (A^3/atom) |
|---|---|---|---|---|---|---|---|---|---|
| TP-LiSbCl4-01 | LiSbCl4 | agm003231612 (NaSbF4) | -54.4 | -73.5 | agm005743014 | -15.0 | agm002157414+agm002205242 | P2_1/c | 24.9 |
| TP-LiSbCl4-02 | LiSbCl4 | agm075662098 (TlAlCl4) | -19.0 | -65.4 | agm005743014 | -7.0 | agm002157414+agm002205242 | P2_1/c | 26.35 |
| TP-LiSbCl4-03 | LiSbCl4 | agm081923399 (CsGaI4) | -8.7 | -64.8 | agm005743014 | -6.4 | agm002157414+agm002205242 | C2/c | 25.74 |
| TP-LiSbCl4-04 | LiSbCl4 | agm085000978 (RbBiF4) | -54.5 | -62.2 | agm005743014 | -3.8 | agm002157414+agm002205242 | P2_1/c | 26.15 |
| TP-LiSbCl4-05 | LiSbCl4 | agm003232459 (CsInI4) | -51.3 | -62.2 | agm005743014 | -3.7 | agm002157414+agm002205242 | P2_1/c | 25.43 |
| TP-LiSbCl4-06 | LiSbCl4 | agm078317478 (CsSbCl4) | -14.1 | -62.1 | agm005743014 | -3.6 | agm002157414+agm002205242 | P1 | 26.27 |
| TP-LiSbCl4-07 | LiSbCl4 | agm078380533 (TlAlCl4) | -3.5 | -62.0 | agm005743014 | -3.5 | agm002157414+agm002205242 | C2/c | 26.06 |
| TP-LiSbCl4-08 | LiSbCl4 | agm077692681 (TlSbF4) | -35.8 | -59.6 | agm005743014 | -1.1 | agm002157414+agm002205242 | P2_1/c | 24.87 |
| TP-LiSbCl4-09 | LiSbCl4 | agm084286613 (NaSbF4) | -48.4 | -58.3 | agm005743014 | +0.2 | agm002157414+agm002205242 | P2_12_12_1 | 25.24 |
| TP-LiSbCl4-10 | LiSbCl4 | agm003232905 (LiGaI4) | -44.9 | -58.1 | agm005743014 | +0.4 | agm002157414+agm002205242 | P2_1/c | 25.61 |
| TP-LiSbCl4-11 | LiSbCl4 | agm077914577 (KBiF4) | -51.0 | -57.6 | agm005743014 | +0.9 | agm002157414+agm002205242 | P2_1 | 25.02 |
| TP-LiSbCl4-12 | LiSbCl4 | agm003222835 (NaYF4) | -47.5 | -57.4 | agm005743014 | +1.1 | agm002157414+agm002205242 | P1 | 25.95 |
| TP-LiSbCl4-13 | LiSbCl4 | agm078792790 (CsSbF4) | -40.0 | -55.9 | agm005743014 | +2.6 | agm002157414+agm002205242 | P2_1/c | 26.48 |
| TP-LiSbCl4-14 | LiSbCl4 | agm020526051 (RbInBr4) | -20.7 | -54.6 | agm005743014 | +3.9 | agm002157414+agm002205242 | P-1 | 26.33 |
| TP-LiSbCl4-15 | LiSbCl4 | agm082661748 (NaYBr4) | -38.5 | -54.0 | agm005743014 | +4.5 | agm002157414+agm002205242 | Pbcn | 25.07 |
| TP-LiSbCl4-16 | LiSbCl4 | agm084497776 (KGaBr4) | -11.3 | -53.6 | agm005743014 | +4.9 | agm002157414+agm002205242 | P2_12_12_1 | 25.8 |
| TP-LiSbCl4-17 | LiSbCl4 | agm005743614 (NaBiCl4) | -54.1 | -51.9 | agm005743014 | +6.6 | agm002157414+agm002205242 | C2/c | 27.67 |
| TP-LiSbCl4-18 | LiSbCl4 | agm083035293 (CsGaCl4) | -50.9 | -51.5 | agm005743014 | +7.0 | agm002157414+agm002205242 | P2_12_12_1 | 27.4 |
| TP-LiSbCl4-19 | LiSbCl4 | agm079169076 (LiGaF4) | -38.7 | -51.3 | agm005743014 | +7.2 | agm002157414+agm002205242 | Pca2_1 | 26.5 |
| TP-LiSbCl4-20 | LiSbCl4 | agm084661913 (NaBiCl4) | -28.9 | -51.3 | agm005743014 | +7.2 | agm002157414+agm002205242 | Fdd2 | 24.61 |
| TP-LiSbCl4-21 | LiSbCl4 | agm085149200 (KBiF4) | -5.8 | -50.6 | agm005743014 | +7.9 | agm002157414+agm002205242 | Pc | 26.61 |
| TP-LiSbCl4-22 | LiSbCl4 | agm002143056 (AgBiCl4) | -24.4 | -50.2 | agm005743014 | +8.3 | agm002157414+agm002205242 | I2_12_12_1 | 25.72 |
| TP-LiSbCl4-23 | LiSbCl4 | agm078016560 (NaGaBr4) | -31.5 | -50.2 | agm005743014 | +8.3 | agm002157414+agm002205242 | P2_12_12_1 | 25.49 |
| TP-LiSbCl4-24 | LiSbCl4 | agm088212697 (RbBiCl4) | -30.1 | -49.9 | agm005743014 | +8.6 | agm002157414+agm002205242 | P2_12_12 | 24.71 |
| TP-LiSbCl4-25 | LiSbCl4 | agm075868128 (RbAlCl4) | -48.9 | -49.3 | agm005743014 | +9.2 | agm002157414+agm002205242 | P2_1/c | 26.92 |
| TP-LiSbCl4-26 | LiSbCl4 | agm082550549 (TlInI4) | -48.4 | -48.8 | agm005743014 | +9.7 | agm002157414+agm002205242 | P2_12_12_1 | 25.87 |
| TP-LiSbCl4-27 | LiSbCl4 | agm003222834 (NaYF4) | -38.6 | -45.8 | agm005743014 | +12.7 | agm002157414+agm002205242 | P3 | 26.24 |
| TP-LiSbCl4-28 | LiSbCl4 | agm077947982 (NaGaF4) | -27.2 | -45.2 | agm005743014 | +13.3 | agm002157414+agm002205242 | Pbcn | 25.03 |
| TP-LiSbCl4-29 | LiSbCl4 | agm013231565 (TlYF4) | -24.8 | -44.6 | agm005743014 | +13.9 | agm002157414+agm002205242 | P1 | 24.97 |
| TP-LiSbCl4-30 | LiSbCl4 | agm003231822 (KGaCl4) | -18.8 | -44.1 | agm005743014 | +14.4 | agm002157414+agm002205242 | P2_1/c | 26.41 |
| TP-LiSbCl4-31 | LiSbCl4 | agm078229809 (TlSbBr4) | -28.9 | -43.9 | agm005743014 | +14.6 | agm002157414+agm002205242 | P2_1/c | 25.87 |
| TP-LiSbCl4-32 | LiSbCl4 | agm083285707 (AgGaBr4) | -37.7 | -43.2 | agm005743014 | +15.3 | agm002157414+agm002205242 | P2_12_12_1 | 27.22 |
| TP-LiSbCl4-33 | LiSbCl4 | agm082993966 (NaBiF4) | -33.4 | -42.5 | agm005743014 | +16.0 | agm002157414+agm002205242 | P-1 | 25.23 |
| TP-LiSbCl4-34 | LiSbCl4 | agm089761403 (RbBiCl4) | -60.3 | -40.2 | agm005743014 | +18.2 | agm002157414+agm002205242 | P2_1/c | 29.95 |
| TP-LiSbCl4-35 | LiSbCl4 | agm085030268 (CsInF4) | -30.2 | -39.8 | agm005743014 | +18.6 | agm002157414+agm002205242 | P-1 | 26.26 |
| TP-LiSbCl4-36 | LiSbCl4 | agm078457792 (KAlCl4) | -29.7 | -38.2 | agm005743014 | +20.3 | agm002157414+agm002205242 | P2_1 | 25.99 |
| TP-NaSbCl4-01 | NaSbCl4 | agm083138585 (CsAlCl4) | -0.0 | -38.0 | agm005743441 | -18.7 | agm003157609+agm002205242 | P-1 | 26.42 |
| TP-LiSbCl4-37 | LiSbCl4 | agm083044527 (RbSbF4) | -43.4 | -37.6 | agm005743014 | +20.9 | agm002157414+agm002205242 | Pna2_1 | 27.12 |
| TP-LiSbCl4-38 | LiSbCl4 | agm082187157 (LiGaF4) | -24.4 | -37.6 | agm005743014 | +20.9 | agm002157414+agm002205242 | Pnma | 25.16 |
| TP-LiSbCl4-39 | LiSbCl4 | agm003231751 (LiInF4) | -8.6 | -37.3 | agm005743014 | +21.2 | agm002157414+agm002205242 | Pbcn | 24.12 |
| TP-LiSbCl4-40 | LiSbCl4 | agm067782374 (CsAlBr4) | -17.4 | -36.8 | agm005743014 | +21.7 | agm002157414+agm002205242 | P1 | 26.18 |
| TP-LiSbCl4-41 | LiSbCl4 | agm003231473 (KAlF4) | -12.4 | -36.7 | agm005743014 | +21.8 | agm002157414+agm002205242 | P2_1 | 26.58 |
| TP-LiSbCl4-42 | LiSbCl4 | agm003232193 (TlSbF4) | -9.6 | -35.9 | agm005743014 | +22.5 | agm002157414+agm002205242 | P2_1/c | 24.23 |
| TP-LiSbCl4-43 | LiSbCl4 | agm088369391 (CsGaBr4) | -33.0 | -35.9 | agm005743014 | +22.6 | agm002157414+agm002205242 | Pna2_1 | 26.27 |
| TP-LiSbCl4-44 | LiSbCl4 | agm003736958 (TlSbCl4) | -18.7 | -35.6 | agm005743014 | +22.9 | agm002157414+agm002205242 | C2/c | 25.45 |
| TP-LiSbCl4-45 | LiSbCl4 | agm078214723 (NaYCl4) | -8.1 | -35.1 | agm005743014 | +23.4 | agm002157414+agm002205242 | P2/c | 24.33 |
| TP-LiSbCl4-46 | LiSbCl4 | agm005733339 (CsSbCl4) | -16.8 | -35.1 | agm005743014 | +23.4 | agm002157414+agm002205242 | C2/m | 26.39 |
| TP-LiSbCl4-47 | LiSbCl4 | agm085182245 (NaInF4) | -31.7 | -35.0 | agm005743014 | +23.4 | agm002157414+agm002205242 | Pc | 25.48 |
| TP-LiSbCl4-48 | LiSbCl4 | agm085286343 (KSbF4) | -43.8 | -34.7 | agm005743014 | +23.8 | agm002157414+agm002205242 | Pc | 28.27 |
| TP-LiSbCl4-49 | LiSbCl4 | agm003212503 (TlBiF4) | -32.3 | -34.5 | agm005743014 | +24.0 | agm002157414+agm002205242 | P2_1 | 26.48 |
| TP-NaSbCl4-02 | NaSbCl4 | agm003231612 (NaSbF4) | -54.1 | -34.3 | agm005743441 | -15.0 | agm003157609+agm002205242 | P2_1/c | 28.2 |
| TP-LiSbCl4-50 | LiSbCl4 | agm003232194 (TlSbF4) | -41.0 | -33.9 | agm005743014 | +24.6 | agm002157414+agm002205242 | P-1 | 27.51 |
| TP-LiBiCl4-01 | LiBiCl4 | agm084160006 (RbAlCl4) | -26.4 | -33.6 | agm003736573 | -13.6 | agm002157414+agm002142653 | P2_1/c | 24.38 |
| TP-LiSbCl4-51 | LiSbCl4 | agm005740896 (AgYCl4) | -21.2 | -33.3 | agm005743014 | +25.2 | agm002157414+agm002205242 | P2/m | 25.64 |
| TP-LiSbCl4-52 | LiSbCl4 | agm003740511 (LiInBr4) | -21.7 | -33.2 | agm005743014 | +25.3 | agm002157414+agm002205242 | C2/m | 25.55 |
| TP-LiSbCl4-53 | LiSbCl4 | agm088446983 (KGaCl4) | -3.6 | -33.1 | agm005743014 | +25.4 | agm002157414+agm002205242 | P-1 | 27.57 |
| TP-LiSbCl4-54 | LiSbCl4 | agm003736714 (NaBiCl4) | -22.4 | -32.4 | agm005743014 | +26.1 | agm002157414+agm002205242 | P2/m | 25.51 |
| TP-NaSbCl4-03 | NaSbCl4 | agm078334508 (CsGaBr4) | -35.0 | -31.6 | agm005743441 | -12.3 | agm003157609+agm002205242 | P2_12_12_1 | 26.51 |
| TP-LiSbCl4-55 | LiSbCl4 | agm005739607 (NaSbBr4) | -23.7 | -31.4 | agm005743014 | +27.1 | agm002157414+agm002205242 | C2 | 25.53 |
| TP-LiBiCl4-02 | LiBiCl4 | agm003231612 (NaSbF4) | -18.4 | -31.2 | agm003736573 | -11.2 | agm002157414+agm002142653 | P2_1/c | 24.71 |
| TP-LiSbCl4-56 | LiSbCl4 | agm087483791 (NaAlCl4) | -34.5 | -31.1 | agm005743014 | +27.3 | agm002157414+agm002205242 | Pnnm | 27.64 |
| TP-LiSbCl4-57 | LiSbCl4 | agm078069920 (CsAlI4) | -39.1 | -30.8 | agm005743014 | +27.7 | agm002157414+agm002205242 | P2_1/c | 28.28 |
| TP-LiSbCl4-58 | LiSbCl4 | agm084475464 (CsGaBr4) | +7.1 | -30.8 | agm005743014 | +27.7 | agm002157414+agm002205242 | P-1 | 26.78 |
| TP-NaSbCl4-04 | NaSbCl4 | agm083044527 (RbSbF4) | -25.2 | -30.6 | agm005743441 | -11.3 | agm003157609+agm002205242 | Pna2_1 | 26.02 |
| TP-NaSbCl4-05 | NaSbCl4 | agm003231822 (KGaCl4) | -37.9 | -30.4 | agm005743441 | -11.1 | agm003157609+agm002205242 | P2_1/c | 26.64 |
| TP-LiSbCl4-59 | LiSbCl4 | agm082822284 (LiGaF4) | +1.9 | -30.3 | agm005743014 | +28.2 | agm002157414+agm002205242 | P2_1 | 25.2 |
| TP-LiSbCl4-60 | LiSbCl4 | agm078239010 (NaInBr4) | -36.0 | -30.3 | agm005743014 | +28.2 | agm002157414+agm002205242 | P-1 | 28.4 |
| TP-LiBiCl4-03 | LiBiCl4 | agm003232459 (CsInI4) | -26.2 | -30.1 | agm003736573 | -10.1 | agm002157414+agm002142653 | P2_1/c | 25.56 |
| TP-LiSbCl4-61 | LiSbCl4 | agm081414147 (CsBiCl4) | -19.0 | -29.9 | agm005743014 | +28.5 | agm002157414+agm002205242 | P2/c | 25.9 |
| TP-LiSbCl4-62 | LiSbCl4 | agm003741334 (AgInF4) | -43.2 | -29.9 | agm005743014 | +28.6 | agm002157414+agm002205242 | P1 | 27.95 |
| TP-LiSbCl4-63 | LiSbCl4 | agm075717463 (RbGaCl4) | -20.6 | -29.9 | agm005743014 | +28.6 | agm002157414+agm002205242 | P2_1/c | 26.03 |
| TP-LiSbCl4-64 | LiSbCl4 | agm003739564 (NaScCl4) | -26.2 | -29.7 | agm005743014 | +28.8 | agm002157414+agm002205242 | C2/c | 26.01 |
| TP-LiSbCl4-65 | LiSbCl4 | agm082986275 (KBiCl4) | -22.8 | -29.2 | agm005743014 | +29.3 | agm002157414+agm002205242 | Pnna | 26.34 |
| TP-LiSbCl4-66 | LiSbCl4 | agm083575773 (TlBiBr4) | -26.9 | -29.1 | agm005743014 | +29.4 | agm002157414+agm002205242 | C2/c | 26.21 |
| TP-LiSbCl4-67 | LiSbCl4 | agm082332867 (TlAlCl4) | -5.5 | -28.7 | agm005743014 | +29.8 | agm002157414+agm002205242 | Pna2_1 | 26.73 |
| TP-LiSbCl4-68 | LiSbCl4 | agm080048502 (KAlCl4) | -33.6 | -28.3 | agm005743014 | +30.2 | agm002157414+agm002205242 | P2_1 | 27.97 |
| TP-LiSbCl4-69 | LiSbCl4 | agm005729801 (KInF4) | -35.6 | -28.2 | agm005743014 | +30.3 | agm002157414+agm002205242 | P1 | 27.38 |
| TP-LiSbCl4-70 | LiSbCl4 | agm077761097 (CsGaI4) | -21.8 | -28.0 | agm005743014 | +30.5 | agm002157414+agm002205242 | P2_1/c | 26.44 |
| TP-LiBiCl4-04 | LiBiCl4 | agm083138585 (CsAlCl4) | -13.2 | -27.8 | agm003736573 | -7.8 | agm002157414+agm002142653 | P-1 | 25.97 |
| TP-LiSbCl4-71 | LiSbCl4 | agm077779806 (KGaCl4) | -17.7 | -27.8 | agm005743014 | +30.7 | agm002157414+agm002205242 | P2_12_12_1 | 27.24 |
| TP-LiSbCl4-72 | LiSbCl4 | agm075849623 (CsGaCl4) | -18.1 | -27.0 | agm005743014 | +31.5 | agm002157414+agm002205242 | Pna2_1 | 27.36 |
| TP-LiSbCl4-73 | LiSbCl4 | agm078953375 (RbGaCl4) | -34.4 | -26.4 | agm005743014 | +32.1 | agm002157414+agm002205242 | Pna2_1 | 28.07 |
| TP-LiSbCl4-74 | LiSbCl4 | agm064032023 (NaInF4) | -15.2 | -26.2 | agm005743014 | +32.2 | agm002157414+agm002205242 | C2/m | 25.42 |
| TP-LiSbCl4-75 | LiSbCl4 | agm082259872 (RbGaCl4) | -27.3 | -26.2 | agm005743014 | +32.3 | agm002157414+agm002205242 | P-1 | 28.9 |
| TP-NaSbCl4-06 | NaSbCl4 | agm077764310 (TlBiF4) | +7.8 | -26.0 | agm005743441 | -6.7 | agm003157609+agm002205242 | P2_1/c | 26.88 |
| TP-LiSbCl4-76 | LiSbCl4 | agm083199798 (NaYF4) | -23.1 | -26.0 | agm005743014 | +32.5 | agm002157414+agm002205242 | P2_1 | 27.34 |
| TP-LiBiCl4-05 | LiBiCl4 | agm003232784 (LiGaBr4) | -16.8 | -25.6 | agm003736573 | -5.6 | agm002157414+agm002142653 | P2_1/c | 25.84 |
| TP-NaSbCl4-07 | NaSbCl4 | agm005729801 (KInF4) | +10.9 | -25.6 | agm005743441 | -6.3 | agm003157609+agm002205242 | P2_1 | 26.81 |
| TP-LiSbCl4-77 | LiSbCl4 | agm003231472 (LiAlCl4) | -27.2 | -25.5 | agm005743014 | +33.0 | agm002157414+agm002205242 | P2_1/c | 26.44 |
| TP-NaSbCl4-08 | NaSbCl4 | agm084286613 (NaSbF4) | -45.9 | -25.4 | agm005743441 | -6.1 | agm003157609+agm002205242 | P2_12_12_1 | 28.03 |
| TP-LiSbCl4-78 | LiSbCl4 | agm078075920 (KAlCl4) | -1.6 | -25.2 | agm005743014 | +33.3 | agm002157414+agm002205242 | P1 | 27.18 |
| TP-LiBiCl4-06 | LiBiCl4 | agm088642654 (TlScCl4) | -8.9 | -24.9 | agm003736573 | -4.9 | agm002157414+agm002142653 | Pbcn | 25.36 |
| TP-LiBiCl4-07 | LiBiCl4 | agm077914577 (KBiF4) | -16.9 | -24.8 | agm003736573 | -4.8 | agm002157414+agm002142653 | P2_1 | 25.31 |
| TP-LiSbCl4-79 | LiSbCl4 | agm083089336 (NaInF4) | -9.9 | -23.9 | agm005743014 | +34.5 | agm002157414+agm002205242 | P2_1/c | 25.61 |
| TP-NaSbCl4-09 | NaSbCl4 | agm075868128 (RbAlCl4) | -20.5 | -23.8 | agm005743441 | -4.5 | agm003157609+agm002205242 | P2_1/c | 27.53 |
| TP-LiSbBr4-01 | LiSbBr4 | agm078380533 (TlAlCl4) | -23.6 | -23.7 | agm005739369 | -19.7 | agm002091979+agm002139075 | C2/c | 29.45 |
| TP-LiSbBr4-02 | LiSbBr4 | agm085286343 (KSbF4) | +11.4 | -23.4 | agm005739369 | -19.4 | agm002091979+agm002139075 | P2/c | 29.24 |
| TP-LiBiCl4-08 | LiBiCl4 | agm082360769 (RbAlI4) | +4.1 | -23.4 | agm003736573 | -3.4 | agm002157414+agm002142653 | P-1 | 24.51 |
| TP-LiBiBr4-01 | LiBiBr4 | agm082912513 (RbGaF4) | -14.9 | -23.4 | agm005739442 | -11.1 | agm002091979+agm002138960 | C2/c | 29.71 |
| TP-LiSbCl4-80 | LiSbCl4 | agm087902309 (KSbF4) | -8.3 | -22.7 | agm005743014 | +35.8 | agm002157414+agm002205242 | Imma | 25.63 |
| TP-NaSbBr4-01 | NaSbBr4 | agm078380533 (TlAlCl4) | -22.1 | -22.0 | agm003733885 | -5.1 | agm003157587+agm002139075 | C2/c | 32.45 |
| TP-LiSbCl4-81 | LiSbCl4 | agm082310455 (CsInI4) | -4.5 | -21.6 | agm005743014 | +36.9 | agm002157414+agm002205242 | P1 | 27.35 |
| TP-LiSbBr4-03 | LiSbBr4 | agm082912513 (RbGaF4) | -16.2 | -21.5 | agm005739369 | -17.4 | agm002091979+agm002139075 | P-1 | 29.37 |
| TP-NaSbCl4-10 | NaSbCl4 | agm077835909 (RbInCl4) | -37.1 | -21.3 | agm005743441 | -2.0 | agm003157609+agm002205242 | P2_1/c | 27.54 |
| TP-NaGaI4-01 | NaGaI4 | agm003232905 (LiGaI4) | +6.1 | -21.3 | agm003736187 | -8.1 | agm003157588+agm003220383 | P2_1/c | 38.38 |
| TP-LiBiCl4-09 | LiBiCl4 | agm088212697 (RbBiCl4) | -12.0 | -20.9 | agm003736573 | -0.9 | agm002157414+agm002142653 | P2_12_12 | 25.17 |
| TP-NaSbCl4-11 | NaSbCl4 | agm085000978 (RbBiF4) | -47.1 | -20.8 | agm005743441 | -1.4 | agm003157609+agm002205242 | P2_1/c | 28.9 |
| TP-LiBiBr4-02 | LiBiBr4 | agm078380533 (TlAlCl4) | -19.7 | -20.6 | agm005739442 | -8.3 | agm002091979+agm002138960 | C2/c | 30.03 |
| TP-LiBiCl4-10 | LiBiCl4 | agm081123431 (TlBiCl4) | -2.8 | -20.6 | agm003736573 | -0.6 | agm002157414+agm002142653 | P2/c | 24.85 |
| TP-NaSbI4-01 | NaSbI4 | agm083138585 (CsAlCl4) | -16.1 | -20.5 | agm005743434 | +6.2 | agm003157588+agm005543490 | P-1 | 38.76 |
| TP-LiSbCl4-82 | LiSbCl4 | agm081964834 (KInBr4) | -12.1 | -20.4 | agm005743014 | +38.1 | agm002157414+agm002205242 | C222_1 | 25.18 |
| TP-LiSbCl4-83 | LiSbCl4 | agm084160006 (RbAlCl4) | -28.0 | -20.2 | agm005743014 | +38.3 | agm002157414+agm002205242 | P2_1/c | 27.56 |
| TP-LiBiCl4-11 | LiBiCl4 | agm003231751 (LiInF4) | -3.4 | -20.2 | agm003736573 | -0.2 | agm002157414+agm002142653 | Pbcn | 24.66 |
| TP-LiSbCl4-84 | LiSbCl4 | agm082492622 (KSbF4) | +9.3 | -20.2 | agm005743014 | +38.3 | agm002157414+agm002205242 | P4_3 | 24.39 |
| TP-NaSbCl4-12 | NaSbCl4 | agm077947982 (NaGaF4) | -39.9 | -20.1 | agm005743441 | -0.8 | agm003157609+agm002205242 | Pbcn | 28.13 |
| TP-NaSbCl4-13 | NaSbCl4 | agm082661748 (NaYBr4) | -39.3 | -20.1 | agm005743441 | -0.7 | agm003157609+agm002205242 | Pbcn | 28.07 |
| TP-LiSbCl4-85 | LiSbCl4 | agm082817793 (LiGaF4) | -4.2 | -19.1 | agm005743014 | +39.4 | agm002157414+agm002205242 | Cmmm | 25.23 |
| TP-LiSbBr4-04 | LiSbBr4 | agm003222835 (NaYF4) | -17.8 | -19.1 | agm005739369 | -15.0 | agm002091979+agm002139075 | P1 | 29.21 |
| TP-NaSbCl4-14 | NaSbCl4 | agm078792790 (CsSbF4) | -19.6 | -18.9 | agm005743441 | +0.4 | agm003157609+agm002205242 | P2_1/c | 26.99 |
| TP-LiSbCl4-86 | LiSbCl4 | agm081906546 (AgGaCl4) | -17.5 | -18.8 | agm005743014 | +39.7 | agm002157414+agm002205242 | C222_1 | 26.03 |
| TP-LiSbCl4-87 | LiSbCl4 | agm084951913 (TlSbF4) | -7.8 | -18.7 | agm005743014 | +39.8 | agm002157414+agm002205242 | C2 | 25.66 |
| TP-NaSbCl4-15 | NaSbCl4 | agm080048502 (KAlCl4) | -28.3 | -18.7 | agm005743441 | +0.6 | agm003157609+agm002205242 | P2_1 | 28.88 |
| TP-NaSbCl4-16 | NaSbCl4 | agm003232784 (LiGaBr4) | -40.8 | -18.6 | agm005743441 | +0.7 | agm003157609+agm002205242 | P2_1/c | 28.11 |
| TP-LiSbCl4-88 | LiSbCl4 | agm074178298 (LiInI4) | -1.4 | -18.6 | agm005743014 | +39.9 | agm002157414+agm002205242 | P1 | 26.04 |
| TP-LiSbCl4-89 | LiSbCl4 | agm073334894 (AgAlCl4) | -13.2 | -18.4 | agm005743014 | +40.1 | agm002157414+agm002205242 | P1 | 27.57 |
| TP-LiBiCl4-12 | LiBiCl4 | agm003222835 (NaYF4) | -20.9 | -18.2 | agm003736573 | +1.8 | agm002157414+agm002142653 | C2 | 26.23 |
| TP-NaGaI4-02 | NaGaI4 | agm082912513 (RbGaF4) | -0.5 | -18.1 | agm003736187 | -5.0 | agm003157588+agm003220383 | P-1 | 38.56 |
| TP-NaSbCl4-17 | NaSbCl4 | agm075662098 (TlAlCl4) | -26.8 | -17.8 | agm005743441 | +1.5 | agm003157609+agm002205242 | P2_1/c | 27.11 |
| TP-NaGaI4-03 | NaGaI4 | agm089761403 (RbBiCl4) | -6.1 | -17.6 | agm003736187 | -4.5 | agm003157588+agm003220383 | P2_1/c | 38.64 |
| TP-NaBiBr4-01 | NaBiBr4 | agm083138585 (CsAlCl4) | -15.0 | -17.3 | agm069058740 | -0.9 | agm003157587+agm002138960 | P-1 | 32.85 |
| TP-NaSbCl4-18 | NaSbCl4 | agm088212697 (RbBiCl4) | -31.7 | -16.8 | agm005743441 | +2.5 | agm003157609+agm002205242 | P2_12_12 | 27.39 |
| TP-LiBiCl4-13 | LiBiCl4 | agm005743614 (NaBiCl4) | -26.8 | -16.5 | agm003736573 | +3.5 | agm002157414+agm002142653 | C2/c | 27.12 |
| TP-LiAlBr4-01 | LiAlBr4 | agm064525236 (TlInCl4) | +3.3 | -16.3 | agm003740512 | -47.9 | agm002091979+agm003220531 | P2_1/m | 28.54 |
| TP-LiBiCl4-14 | LiBiCl4 | agm082187157 (LiGaF4) | -8.7 | -14.7 | agm003736573 | +5.3 | agm002157414+agm002142653 | Pnma | 25.49 |
| TP-LiBiBr4-03 | LiBiBr4 | agm003222835 (NaYF4) | -12.5 | -14.5 | agm005739442 | -2.2 | agm002091979+agm002138960 | P1 | 29.81 |
| TP-LiBiCl4-15 | LiBiCl4 | agm085000978 (RbBiF4) | -12.6 | -14.5 | agm003736573 | +5.5 | agm002157414+agm002142653 | P2_1/c | 25.88 |
| TP-LiBiCl4-16 | LiBiCl4 | agm003211813 (LiBiF4) | +14.8 | -14.3 | agm003736573 | +5.7 | agm002157414+agm002142653 | I4_1/a | 23.58 |
| TP-NaSbCl4-19 | NaSbCl4 | agm003231472 (LiAlCl4) | -3.3 | -14.2 | agm005743441 | +5.1 | agm003157609+agm002205242 | P2_1/c | 29.09 |
| TP-NaSbBr4-02 | NaSbBr4 | agm088642654 (TlScCl4) | -5.4 | -14.1 | agm003733885 | +2.8 | agm003157587+agm002139075 | Pbcn | 32.11 |
| TP-LiSbCl4-90 | LiSbCl4 | agm003232962 (NaAlCl4) | -2.2 | -14.0 | agm005743014 | +44.5 | agm002157414+agm002205242 | P2_12_12_1 | 25.62 |
| TP-NaInBr4-01 | NaInBr4 | agm005739607 (NaSbBr4) | +8.8 | -13.9 | agm003229838 | +18.6 | agm003157587+agm002201463 | C2/m | 29.82 |
| TP-NaInBr4-02 | NaInBr4 | agm005742973 (LiScCl4) | +11.8 | -13.9 | agm003229838 | +18.6 | agm003157587+agm002201463 | P2/m | 29.7 |
| TP-LiAlBr4-02 | LiAlBr4 | agm005733339 (CsSbCl4) | +12.9 | -13.7 | agm003740512 | -45.4 | agm002091979+agm003220531 | P-1 | 28.8 |
| TP-LiScCl4-01 | LiScCl4 | agm089761403 (RbBiCl4) | -4.7 | -13.7 | agm005742973 | -3.5 | agm002157414+agm002142959 | P2_1/c | 22.32 |
| TP-NaSbCl4-20 | NaSbCl4 | agm003231618 (RbGaBr4) | -9.4 | -13.7 | agm005743441 | +5.6 | agm003157609+agm002205242 | Pna2_1 | 26.95 |
| TP-NaGaI4-04 | NaGaI4 | agm064525236 (TlInCl4) | -2.2 | -13.6 | agm003736187 | -0.4 | agm003157588+agm003220383 | P2_1/m | 38.9 |
| TP-NaSbBr4-03 | NaSbBr4 | agm089761403 (RbBiCl4) | -8.4 | -13.6 | agm003733885 | +3.3 | agm003157587+agm002139075 | P2_1/c | 31.77 |
| TP-LiBiCl4-17 | LiBiCl4 | agm078457792 (KAlCl4) | -9.1 | -13.6 | agm003736573 | +6.5 | agm002157414+agm002142653 | P2_1 | 26.07 |
| TP-NaSbCl4-21 | NaSbCl4 | agm085030268 (CsInF4) | -23.9 | -13.5 | agm005743441 | +5.8 | agm003157609+agm002205242 | P-1 | 28.27 |
| TP-LiBiCl4-18 | LiBiCl4 | agm002143056 (AgBiCl4) | -12.8 | -13.4 | agm003736573 | +6.5 | agm002157414+agm002142653 | Imma | 25.87 |
| TP-NaSbCl4-22 | NaSbCl4 | agm003222836 (NaYF4) | -24.9 | -13.0 | agm005743441 | +6.3 | agm003157609+agm002205242 | P1 | 27.5 |
| TP-LiSbCl4-91 | LiSbCl4 | agm082246814 (KBiF4) | -22.0 | -12.5 | agm005743014 | +46.0 | agm002157414+agm002205242 | P2_12_12_1 | 27.57 |
| TP-LiBiCl4-19 | LiBiCl4 | agm077947982 (NaGaF4) | -5.6 | -12.4 | agm003736573 | +7.6 | agm002157414+agm002142653 | Pbcn | 25.27 |
| TP-NaGaI4-05 | NaGaI4 | agm084160006 (RbAlCl4) | +0.3 | -12.4 | agm003736187 | +0.8 | agm003157588+agm003220383 | P2_1/c | 39.11 |
| TP-LiSbCl4-92 | LiSbCl4 | agm005755324 (TlGaCl4) | +8.2 | -12.2 | agm005743014 | +46.3 | agm002157414+agm002205242 | P2_1/c | 26.51 |
| TP-LiBiCl4-20 | LiBiCl4 | agm082492622 (KSbF4) | +13.8 | -12.1 | agm003736573 | +7.9 | agm002157414+agm002142653 | P4_3 | 23.76 |
| TP-LiSbCl4-93 | LiSbCl4 | agm088359334 (CsInF4) | -21.5 | -12.0 | agm005743014 | +46.5 | agm002157414+agm002205242 | P1 | 30.9 |
| TP-NaBiBr4-02 | NaBiBr4 | agm082661748 (NaYBr4) | -11.8 | -11.9 | agm069058740 | +4.4 | agm003157587+agm002138960 | Pbcn | 32.78 |
| TP-LiBiCl4-21 | LiBiCl4 | agm082086605 (KGaI4) | -11.9 | -11.9 | agm003736573 | +8.1 | agm002157414+agm002142653 | P2_12_12_1 | 26.05 |
| TP-NaSbBr4-04 | NaSbBr4 | agm003232459 (CsInI4) | -3.8 | -11.7 | agm003733885 | +5.1 | agm003157587+agm002139075 | P2_1/c | 32.08 |
| TP-NaInBr4-03 | NaInBr4 | agm088642654 (TlScCl4) | +12.2 | -11.7 | agm003229838 | +20.8 | agm003157587+agm002201463 | Pbcn | 30.05 |
| TP-LiScCl4-02 | LiScCl4 | agm083138585 (CsAlCl4) | -15.8 | -11.0 | agm005742973 | -0.8 | agm002157414+agm002142959 | P-1 | 22.87 |
| TP-LiBiCl4-22 | LiBiCl4 | agm083044527 (RbSbF4) | -18.6 | -11.0 | agm003736573 | +9.0 | agm002157414+agm002142653 | Pna2_1 | 26.68 |
| TP-LiBiCl4-23 | LiBiCl4 | agm077692681 (TlSbF4) | +11.2 | -10.7 | agm003736573 | +9.3 | agm002157414+agm002142653 | P2_1/c | 23.59 |
| TP-NaSbCl4-23 | NaSbCl4 | agm077692681 (TlSbF4) | -30.6 | -10.7 | agm005743441 | +8.6 | agm003157609+agm002205242 | P2_1/c | 27.77 |
| TP-NaInBr4-04 | NaInBr4 | agm003736714 (NaBiCl4) | +8.8 | -10.7 | agm003229838 | +21.8 | agm003157587+agm002201463 | P2/m | 30.05 |
| TP-LiAlBr4-03 | LiAlBr4 | agm082993966 (NaBiF4) | +1.6 | -10.5 | agm003740512 | -42.1 | agm002091979+agm003220531 | P-1 | 29.38 |
| TP-LiBiBr4-04 | LiBiBr4 | agm003232459 (CsInI4) | -5.0 | -10.4 | agm005739442 | +1.9 | agm002091979+agm002138960 | P2_1/c | 29.84 |
| TP-NaSbCl4-24 | NaSbCl4 | agm078016560 (NaGaBr4) | -32.0 | -10.4 | agm005743441 | +8.9 | agm003157609+agm002205242 | P2_12_12_1 | 28.11 |
| TP-NaSbCl4-25 | NaSbCl4 | agm082492622 (KSbF4) | -18.1 | -10.3 | agm005743441 | +9.0 | agm003157609+agm002205242 | P4_3 | 26.63 |
| TP-NaSbCl4-26 | NaSbCl4 | agm005755324 (TlGaCl4) | -11.4 | -10.2 | agm005743441 | +9.1 | agm003157609+agm002205242 | P-1 | 26.64 |
| TP-NaGaI4-06 | NaGaI4 | agm003232962 (NaAlCl4) | +0.9 | -10.1 | agm003736187 | +3.1 | agm003157588+agm003220383 | P2_12_12_1 | 37.7 |
| TP-LiSbCl4-94 | LiSbCl4 | agm081555238 (CsSbCl4) | -31.4 | -10.1 | agm005743014 | +48.4 | agm002157414+agm002205242 | P2_1/c | 29.64 |
| TP-LiBiBr4-05 | LiBiBr4 | agm082661748 (NaYBr4) | -2.8 | -9.9 | agm005739442 | +2.4 | agm002091979+agm002138960 | Pbcn | 29.58 |
| TP-LiYBr4-01 | LiYBr4 | agm013231565 (TlYF4) | +13.9 | -9.7 | agm003742437 | -10.6 | agm002091979+agm002139115 | Pnma | 26.78 |
| TP-NaSbCl4-27 | NaSbCl4 | agm078214723 (NaYCl4) | -22.4 | -9.7 | agm005743441 | +9.7 | agm003157609+agm002205242 | P2/c | 27.37 |
| TP-NaBiI4-01 | NaBiI4 | agm078380533 (TlAlCl4) | -8.7 | -9.6 | agm003740643 | +5.7 | agm003157588+agm002154521 | C2/c | 39.83 |
| TP-NaSbCl4-28 | NaSbCl4 | agm077761097 (CsGaI4) | -4.2 | -9.6 | agm005743441 | +9.7 | agm003157609+agm002205242 | P2_1/c | 27.42 |
| TP-NaBiBr4-03 | NaBiBr4 | agm003222835 (NaYF4) | -11.2 | -9.4 | agm069058740 | +6.9 | agm003157587+agm002138960 | P1 | 32.54 |
| TP-NaSbBr4-05 | NaSbBr4 | agm003232784 (LiGaBr4) | -4.8 | -9.3 | agm003733885 | +7.6 | agm003157587+agm002139075 | P2_1/c | 32.29 |
| TP-NaSbCl4-29 | NaSbCl4 | agm003232194 (TlSbF4) | -13.6 | -9.0 | agm005743441 | +10.3 | agm003157609+agm002205242 | P-1 | 28.21 |
| TP-LiSbCl4-95 | LiSbCl4 | agm002218760 (LiAlF4) | -6.4 | -9.0 | agm005743014 | +49.5 | agm002157414+agm002205242 | P1 | 32.01 |
| TP-NaSbI4-02 | NaSbI4 | agm089761403 (RbBiCl4) | +4.4 | -8.9 | agm005743434 | +17.7 | agm003157588+agm005543490 | P2_1/c | 38.05 |
| TP-NaScBr4-01 | NaScBr4 | agm082661748 (NaYBr4) | -12.9 | -8.7 | agm005739585 | +12.2 | agm003157587+agm002184236 | Pbcn | 29.68 |
| TP-LiSbBr4-05 | LiSbBr4 | agm082661748 (NaYBr4) | -4.1 | -8.6 | agm005739369 | -4.6 | agm002091979+agm002139075 | Pbcn | 28.87 |
| TP-NaSbCl4-30 | NaSbCl4 | agm005743614 (NaBiCl4) | -48.2 | -8.6 | agm005743441 | +10.7 | agm003157609+agm002205242 | C2/c | 30.55 |
| TP-NaGaI4-07 | NaGaI4 | agm078792790 (CsSbF4) | +10.9 | -8.5 | agm003736187 | +4.7 | agm003157588+agm003220383 | P2_1/c | 39.53 |
| TP-NaSbI4-03 | NaSbI4 | agm088642654 (TlScCl4) | +6.1 | -8.4 | agm005743434 | +18.2 | agm003157588+agm005543490 | Pbcn | 39.12 |
| TP-LiSbBr4-06 | LiSbBr4 | agm003232459 (CsInI4) | -5.7 | -8.2 | agm005739369 | -4.2 | agm002091979+agm002139075 | P2_1/c | 29.16 |
| TP-NaBiCl4-01 | NaBiCl4 | agm082086605 (KGaI4) | +7.3 | -8.1 | agm088474086 | -16.4 | agm003157609+agm002142653 | P2_12_12_1 | 26.11 |
| TP-LiBiBr4-06 | LiBiBr4 | agm003232784 (LiGaBr4) | +0.7 | -8.0 | agm005739442 | +4.3 | agm002091979+agm002138960 | P2_1/c | 29.86 |
| TP-NaBiBr4-04 | NaBiBr4 | agm003232459 (CsInI4) | -3.1 | -8.0 | agm069058740 | +8.4 | agm003157587+agm002138960 | P2_1/c | 32.69 |
| TP-LiBiCl4-24 | LiBiCl4 | agm083089336 (NaInF4) | +10.8 | -7.9 | agm003736573 | +12.1 | agm002157414+agm002142653 | P2_1/c | 24.87 |
| TP-LiSbCl4-96 | LiSbCl4 | agm083289675 (TlAlCl4) | -6.7 | -7.8 | agm005743014 | +50.7 | agm002157414+agm002205242 | P2/c | 26.45 |
| TP-LiInBr4-01 | LiInBr4 | agm082912513 (RbGaF4) | +11.7 | -7.5 | agm003740511 | -11.4 | agm002091979+agm002201463 | P-1 | 26.97 |
| TP-NaSbI4-04 | NaSbI4 | agm003740643 (NaBiI4) | -2.9 | -7.4 | agm005743434 | +19.2 | agm003157588+agm005543490 | P2/m | 38.16 |
| TP-LiBiCl4-25 | LiBiCl4 | agm084497776 (KGaBr4) | +0.8 | -7.2 | agm003736573 | +12.8 | agm002157414+agm002142653 | P2_12_12_1 | 25.19 |
| TP-LiInCl4-01 | LiInCl4 | agm085030268 (CsInF4) | -4.4 | -7.0 | agm003740637 | -10.0 | agm002157414+agm002142775 | P-1 | 23.21 |
| TP-LiYBr4-02 | LiYBr4 | agm077835909 (RbInCl4) | -4.7 | -7.0 | agm003742437 | -7.9 | agm002091979+agm002139115 | C2/c | 28.9 |
| TP-LiScBr4-01 | LiScBr4 | agm082912513 (RbGaF4) | -6.6 | -6.8 | agm005739333 | -10.9 | agm002091979+agm002184236 | P-1 | 26.47 |
| TP-NaSbCl4-31 | NaSbCl4 | agm003231833 (RbAlF4) | +9.9 | -6.6 | agm005743441 | +12.7 | agm003157609+agm002205242 | Pmc2_1 | 28.69 |
| TP-LiScCl4-03 | LiScCl4 | agm081414147 (CsBiCl4) | +2.5 | -6.5 | agm005742973 | +3.7 | agm002157414+agm002142959 | P2/c | 21.89 |
| TP-NaBiBr4-05 | NaBiBr4 | agm003232784 (LiGaBr4) | -6.9 | -6.5 | agm069058740 | +9.9 | agm003157587+agm002138960 | P2_1/c | 32.89 |
| TP-LiAlBr4-04 | LiAlBr4 | agm083285707 (AgGaBr4) | -2.5 | -6.3 | agm003740512 | -37.9 | agm002091979+agm003220531 | P2_12_12_1 | 29.73 |
| TP-LiBiCl4-26 | LiBiCl4 | agm079169076 (LiGaF4) | -2.0 | -6.2 | agm003736573 | +13.8 | agm002157414+agm002142653 | Pbcn | 25.43 |
| TP-LiSbBr4-07 | LiSbBr4 | agm003232784 (LiGaBr4) | -3.0 | -6.2 | agm005739369 | -2.2 | agm002091979+agm002139075 | P2_1/c | 29.25 |
| TP-NaInCl4-01 | NaInCl4 | agm082661748 (NaYBr4) | -6.1 | -6.2 | agm069515018 | +1.4 | agm003157609+agm002142775 | Pbcn | 26.21 |
| TP-NaSbCl4-32 | NaSbCl4 | agm085149200 (KBiF4) | -13.9 | -6.1 | agm005743441 | +13.2 | agm003157609+agm002205242 | Pc | 26.99 |
| TP-NaSbBr4-06 | NaSbBr4 | agm003231822 (KGaCl4) | -4.4 | -6.1 | agm003733885 | +10.8 | agm003157587+agm002139075 | P2_1/c | 32.35 |
| TP-NaInCl4-02 | NaInCl4 | agm077947982 (NaGaF4) | -2.0 | -6.0 | agm069515018 | +1.6 | agm003157609+agm002142775 | Pbcn | 25.87 |
| TP-NaSbCl4-33 | NaSbCl4 | agm085182245 (NaInF4) | -30.1 | -6.0 | agm005743441 | +13.3 | agm003157609+agm002205242 | Pnc2 | 28.45 |
| TP-NaSbI4-05 | NaSbI4 | agm005739607 (NaSbBr4) | -1.2 | -5.9 | agm005743434 | +20.7 | agm003157588+agm005543490 | C2/m | 38.92 |
| TP-LiSbBr4-08 | LiSbBr4 | agm078792790 (CsSbF4) | +9.3 | -5.9 | agm005739369 | -1.9 | agm002091979+agm002139075 | P2_1/c | 29.53 |
| TP-NaSbCl4-34 | NaSbCl4 | agm083089336 (NaInF4) | -21.6 | -5.8 | agm005743441 | +13.5 | agm003157609+agm002205242 | P2_1/c | 27.6 |
| TP-NaGaI4-08 | NaGaI4 | agm078953375 (RbGaCl4) | +4.4 | -5.8 | agm003736187 | +7.4 | agm003157588+agm003220383 | Pna2_1 | 39.9 |
| TP-NaSbCl4-35 | NaSbCl4 | agm082187157 (LiGaF4) | -25.8 | -5.8 | agm005743441 | +13.5 | agm003157609+agm002205242 | Pnma | 28.04 |
| TP-NaBiBr4-06 | NaBiBr4 | agm089761403 (RbBiCl4) | -2.0 | -5.7 | agm069058740 | +10.7 | agm003157587+agm002138960 | P2_1/c | 32.29 |
| TP-NaSbCl4-36 | NaSbCl4 | agm003222835 (NaYF4) | -3.1 | -5.6 | agm005743441 | +13.7 | agm003157609+agm002205242 | P1 | 25.15 |
| TP-NaBiCl4-02 | NaBiCl4 | agm077914577 (KBiF4) | +4.2 | -5.5 | agm088474086 | -13.8 | agm003157609+agm002142653 | P2_1 | 28.14 |
| TP-LiBiCl4-27 | LiBiCl4 | agm005742973 (LiScCl4) | -1.2 | -5.1 | agm003736573 | +14.9 | agm002157414+agm002142653 | P2/m | 26.1 |
| TP-LiInCl4-02 | LiInCl4 | agm082912513 (RbGaF4) | -2.2 | -5.1 | agm003740637 | -8.1 | agm002157414+agm002142775 | P-1 | 23.55 |
| TP-NaInBr4-05 | NaInBr4 | agm002143056 (AgBiCl4) | +7.5 | -5.1 | agm003229838 | +27.4 | agm003157587+agm002201463 | Imma | 31.13 |
| TP-NaSbCl4-37 | NaSbCl4 | agm003232193 (TlSbF4) | +8.6 | -5.1 | agm005743441 | +14.2 | agm003157609+agm002205242 | P2_1/c | 25.17 |
| TP-LiSbCl4-97 | LiSbCl4 | agm077764310 (TlBiF4) | -25.6 | -4.8 | agm005743014 | +53.6 | agm002157414+agm002205242 | P2_1/c | 29.91 |
| TP-LiBiCl4-28 | LiBiCl4 | agm003736714 (NaBiCl4) | -7.3 | -4.8 | agm003736573 | +15.2 | agm002157414+agm002142653 | P2/m | 25.88 |
| TP-NaInBr4-06 | NaInBr4 | agm005743614 (NaBiCl4) | +5.7 | -4.8 | agm003229838 | +27.7 | agm003157587+agm002201463 | C2/c | 31.12 |
| TP-LiSbCl4-98 | LiSbCl4 | agm077835909 (RbInCl4) | -10.1 | -4.8 | agm005743014 | +53.7 | agm002157414+agm002205242 | P2_1/c | 27.45 |
| TP-NaSbI4-06 | NaSbI4 | agm085030268 (CsInF4) | -0.6 | -4.8 | agm005743434 | +21.9 | agm003157588+agm005543490 | P-1 | 39.51 |
| TP-LiBiBr4-07 | LiBiBr4 | agm002143056 (AgBiCl4) | -3.7 | -4.7 | agm005739442 | +7.6 | agm002091979+agm002138960 | Imma | 30.22 |
| TP-NaSbCl4-38 | NaSbCl4 | agm082822284 (LiGaF4) | -5.4 | -4.7 | agm005743441 | +14.6 | agm003157609+agm002205242 | Pna2_1 | 27.88 |
| TP-LiScCl4-04 | LiScCl4 | agm077947982 (NaGaF4) | +5.7 | -4.7 | agm005742973 | +5.6 | agm002157414+agm002142959 | Pbcn | 22.04 |
| TP-NaSbCl4-39 | NaSbCl4 | agm078069920 (CsAlI4) | -32.1 | -4.6 | agm005743441 | +14.7 | agm003157609+agm002205242 | P2_1/c | 29.46 |
| TP-LiBiCl4-29 | LiBiCl4 | agm003740511 (LiInBr4) | -5.3 | -4.5 | agm003736573 | +15.5 | agm002157414+agm002142653 | C2/m | 25.99 |
| TP-NaSbBr4-07 | NaSbBr4 | agm077947982 (NaGaF4) | -1.6 | -4.5 | agm003733885 | +12.4 | agm003157587+agm002139075 | Pbcn | 32.42 |
| TP-LiAlBr4-05 | LiAlBr4 | agm089325837 (RbAlCl4) | +2.5 | -4.4 | agm003740512 | -36.1 | agm002091979+agm003220531 | Pmn2_1 | 29.75 |
| TP-NaSbCl4-40 | NaSbCl4 | agm003232459 (CsInI4) | -32.4 | -4.4 | agm005743441 | +14.9 | agm003157609+agm002205242 | P2_1/c | 28.95 |
| TP-NaGaI4-09 | NaGaI4 | agm083138585 (CsAlCl4) | -3.2 | -4.3 | agm003736187 | +8.9 | agm003157588+agm003220383 | P-1 | 38.69 |
| TP-LiBiCl4-30 | LiBiCl4 | agm078229809 (TlSbBr4) | -3.7 | -4.3 | agm003736573 | +15.7 | agm002157414+agm002142653 | P2_1/c | 26.02 |
| TP-LiBiCl4-31 | LiBiCl4 | agm069948232 (RbBiCl4) | -4.2 | -4.3 | agm003736573 | +15.7 | agm002157414+agm002142653 | C2/c | 25.84 |
| TP-LiBiCl4-32 | LiBiCl4 | agm082912513 (RbGaF4) | -12.7 | -4.2 | agm003736573 | +15.8 | agm002157414+agm002142653 | P-1 | 26.8 |
| TP-NaSbCl4-41 | NaSbCl4 | agm082912513 (RbGaF4) | -13.1 | -4.2 | agm005743441 | +15.1 | agm003157609+agm002205242 | P-1 | 28.26 |
| TP-NaSbBr4-08 | NaSbBr4 | agm021417380 (CsInI4) | +0.6 | -4.1 | agm003733885 | +12.8 | agm003157587+agm002139075 | P-1 | 32.36 |
| TP-NaGaI4-10 | NaGaI4 | agm005729801 (KInF4) | +14.3 | -3.8 | agm003736187 | +9.3 | agm003157588+agm003220383 | P1 | 38.93 |
| TP-NaScBr4-02 | NaScBr4 | agm077947982 (NaGaF4) | -4.0 | -3.8 | agm005739585 | +17.1 | agm003157587+agm002184236 | Pbcn | 29.39 |
| TP-NaBiCl4-03 | NaBiCl4 | agm077947982 (NaGaF4) | +3.2 | -3.6 | agm088474086 | -12.0 | agm003157609+agm002142653 | Pbcn | 28.48 |
| TP-NaInBr4-07 | NaInBr4 | agm083575773 (TlBiBr4) | +5.9 | -3.6 | agm003229838 | +28.9 | agm003157587+agm002201463 | C2/c | 31.5 |
| TP-LiScCl4-05 | LiScCl4 | agm005743614 (NaBiCl4) | -5.8 | -3.6 | agm005742973 | +6.7 | agm002157414+agm002142959 | C2/c | 22.77 |
| TP-NaSbCl4-42 | NaSbCl4 | agm003231751 (LiInF4) | -12.5 | -3.5 | agm005743441 | +15.8 | agm003157609+agm002205242 | Pbcn | 26.92 |
| TP-LiInCl4-03 | LiInCl4 | agm082661748 (NaYBr4) | +4.7 | -3.5 | agm003740637 | -6.5 | agm002157414+agm002142775 | Pbcn | 23.11 |
| TP-NaScBr4-03 | NaScBr4 | agm089761403 (RbBiCl4) | -1.1 | -3.4 | agm005739585 | +17.5 | agm003157587+agm002184236 | P2_1/c | 28.79 |
| TP-NaSbBr4-09 | NaSbBr4 | agm078229809 (TlSbBr4) | -2.5 | -3.4 | agm003733885 | +13.5 | agm003157587+agm002139075 | P2_1/c | 32.6 |
| TP-LiScCl4-06 | LiScCl4 | agm005739607 (NaSbBr4) | -2.2 | -3.4 | agm005742973 | +6.9 | agm002157414+agm002142959 | C2/m | 22.57 |
| TP-LiSbCl4-99 | LiSbCl4 | agm088618043 (RbSbF4) | -12.1 | -3.3 | agm005743014 | +55.1 | agm002157414+agm002205242 | Pna2_1 | 28.15 |
| TP-NaSbCl4-43 | NaSbCl4 | agm077798528 (RbInI4) | +11.9 | -3.2 | agm005743441 | +16.1 | agm003157609+agm002205242 | P2_1 | 29.54 |
| TP-NaSbI4-07 | NaSbI4 | agm003212503 (TlBiF4) | -0.5 | -3.2 | agm005743434 | +23.4 | agm003157588+agm005543490 | P2_1/m | 38.96 |
| TP-LiBiCl4-33 | LiBiCl4 | agm005739607 (NaSbBr4) | -4.5 | -3.0 | agm003736573 | +17.0 | agm002157414+agm002142653 | C2/m | 25.43 |

Survivors of the identity check that are NOT hits (dE > -3 meV/atom or unconverged): 902; lowest ten:

| composition | prototype source | dE MACE+D3 | converged |
|---|---|---|---|
| LiAlBr4 | agm003231822 (KGaCl4) | -10.7 | False |
| LiSbCl4 | agm003232297 (KScF4) | -10.0 | False |
| LiSbCl4 | agm079171113 (TlGaBr4) | -7.2 | False |
| NaBiCl4 | agm003232459 (CsInI4) | -3.0 | True |
| NaBiCl4 | agm003231612 (NaSbF4) | -2.9 | True |
| LiAlCl4 | agm003232962 (NaAlCl4) | -2.8 | True |
| NaBiCl4 | agm082661748 (NaYBr4) | -2.8 | True |
| NaBiCl4 | agm088642654 (TlScCl4) | -2.8 | True |
| NaSbCl4 | agm020526051 (RbInBr4) | -2.7 | True |
| LiAlBr4 | agm076172157 (AgGaI4) | -2.7 | True |

## Sanity anchors

- (a) NaBiBr4 on the NaBiCl4/KBiBr4 prototype (pool representative agm005743614) vs agm069058740, MACE+D3: **+2.6 meV/atom** (skeptic cycle 2: +2.6). Row status: survives_identity; cheap no-D3 dE -8.28 (skeptic: -8.1 at fmax 0.01).
- (b) NaBiCl4 agm088474086 vs agm005743614, MACE+D3: **-12.6 meV/atom** (skeptic cycle 2: -12.6), i.e. agm088474086 is lower.

## Wall time

- This run: 392.6 min; all runs of this script so far (resumed runs skip finished work): 436.2 min over 4 run(s). 8 single-thread workers.

## Caveats

- In-sample ML potential: MACE-MPA-0 was trained on Alexandria/MPtrj, so the existing entries are in its training set and the transplants are not. A margin of a few meV/atom between an in-sample and an out-of-sample structure is within model error.
- D3 damping: Becke-Johnson with PBE parameters (mace_mp default, torch-dftd). MACE was not trained with D3; D3 shifts relative energies of these halides by 10-20 meV/atom, and zero-damping D3 (as in the existing QE inputs) was not tested.
- Alexandria PBE hull distances are not comparable to MACE or MACE+D3 energies. Decomposition energies here use only AX + MX3 (single lowest-hull polymorph each); other ternaries (e.g. A3MX6) are not considered, so E_decomp is not a hull distance.
- No phonons, no finite temperature, no symmetry constraints; FIRE can stall on saddle points. 'Converged' means fmax < 0.01 eV/A within 600 steps.
- Triage is done without D3 (15 meV/atom window, fmax 0.05, <= 150 steps): a prototype that only becomes competitive once dispersion is on can be missed, and only the three lowest no-D3 existing entries are re-relaxed with D3 unless the chemistry has a hit.
- A weak baseline makes a large dE cheap: check the Alexandria hull distance of the best existing entry (e.g. NaSbCl4 has only two listed entries, at +46 and +148 meV/atom PBE). A transplant can also relax far from its source prototype; the CIF, not the prototype name, is the proposed model.
- Cells limited to <= 24 atoms; 'not listed' means not listed in Alexandria PBE for that reduced formula. MP, OQMD, GNoME, ICSD and literature were not checked.
- Per-row data: `discovery/results/transplant_cycle3.csv`; raw relaxations: `discovery/results/transplant_cycle3/{cheap,d3}/`.
