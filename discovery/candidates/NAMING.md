# Candidate naming: read before using any TP-NaBiCl4 folder

The full overnight scan (2026-10-13 02:07) renumbered candidates and overwrote folders written by the first partial scan.

| Name used in EVIDENCE / SCOREBOARD / overnight job files | Parent prototype | Full-scan folder now holding it |
|---|---|---|
| "TP-NaBiCl4-01" (DFT-tested, lowest so far, 12 atoms) | KBiF4-type, agm077914577 | `TP-NaBiCl4-02/` |
| "TPscan01" rival (24 atoms, NOT yet DFT-tested) | KGaI4-type, agm082086605 | `TP-NaBiCl4-01/` |
| "TPscan03" rival (NOT yet DFT-tested) | NaGaF4-type, agm077947982 | `TP-NaBiCl4-03/` |

The DFT starting geometry actually used for the tested model is `discovery/results/skeptic3/relax/TP-NaBiCl4-01_d3.cif`; its relaxed result is `overnight/night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.out`. Always cite the parent prototype ID, not the folder number.
