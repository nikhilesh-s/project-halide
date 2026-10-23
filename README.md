# Project Halide: code for the predicted NaSbBr4 and LiSbBr4 structure models

This archive holds the code behind the predicted structure models TP-NaSbBr4-01 (NaSbBr4) and TP-LiSbBr4-01 (LiSbBr4). Neither compound has a reported crystal structure. Both structures are predicted, proposed models; none has been synthesized or measured, and no electrolyte is claimed.

A machine-learned potential (MACE-MPA-0) screened structure transplants. Quantum ESPRESSO PBE-D3 then checked each hit against the listed Alexandria structures and the competing phases. Seven compositions pass the energy gates. Two models (NaSbBr4 and LiSbBr4) were taken through convergence, dispersion and phonon checks.

## How it maps to the paper and data

- Raw DFT inputs and outputs: NOMAD dataset "Predicted crystal structure models for NaSbBr4 and LiSbBr4 from a prototype-transplant search (PBE-D3), with controls and negative results (v2)" (v2), doi:10.17172/nomad.j0ce-ctz1 (project page: https://nomad-lab.eu/prod/v1/gui/v2/projects/Wozh3Uq0Rb-XVLm9aW3XLQ). They are not in this archive.
- Earlier validation dataset (v1, not a discovery claim): 10.17172/nomad.7c09-m7v2.
- Paper numbers: `paper/claims.py` reads `discovery/results/verdict.json` and the restored DFT outputs and writes `paper/claims.json`. `paper/fill_preprint.py` fills the placeholders in `paper/preprint_v2.md`. `paper/figures/make_figures.py` redraws the figures and tables.
- Steps to rerun everything: `REPRODUCE.md`.

## Gates

- G1: the ML hit survives a rattle-and-relax check (MACE + D3).
- G2: DFT puts it below every listed Alexandria structure of the same composition (PBE-D3, same settings).
- G3: DFT puts it below all tested competing phases (binaries and listed ternaries).
- G4: the G2 and G3 signs hold at a higher cutoff or denser k-mesh.
- G5: DFT phonons show no imaginary modes.

`discovery/scripts/verdict.py` applies G2 to G4 to the pw.x outputs.

## Layout

- `discovery/scripts/`: screening (`transplant_scan.py`, `screen_mace.py`), gate logic (`verdict.py`), stability window (`window.py`).
- `discovery/candidates/`, `discovery/results/`: candidate structures, ML screening results, and the gate and window reports.
- `overnight/`: DFT job generators (`make_*.py`), the batch runner (`run_night.py`), the QE build script (`setup.sh`) and pseudopotential checksums (`pseudo_md5.txt`).
- `anish_packet/`: the cluster workflow for phonons and convergence checks (`go.sh`, `PLAN`, per-stage scripts).
- `paper/`: claim extraction, preprint template and figure scripts.

## License and citation

Code: MIT (`LICENSE`). Cite this archive with `CITATION.cff`, and cite the NOMAD dataset for the data.
