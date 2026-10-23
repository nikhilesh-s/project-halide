# Reproduce

This archive holds code and small reports. It does not hold the raw DFT outputs or the Quantum ESPRESSO pseudopotentials. The raw outputs are in the NOMAD dataset "Predicted crystal structure models for NaSbBr4 and LiSbBr4 from a prototype-transplant search (PBE-D3), with controls and negative results (v2)", doi:10.17172/nomad.j0ce-ctz1 (project page: https://nomad-lab.eu/prod/v1/gui/v2/projects/Wozh3Uq0Rb-XVLm9aW3XLQ). Validation dataset v1: 10.17172/nomad.7c09-m7v2. Run every command from the archive root.

1. Environment. Use Python 3.11 and install the pinned packages: `pip install -r requirements-ml.lock phonopy==4.8.3`. `discovery/ENV.md` lists the exact versions used (pymatgen 2026.9.24, ase 3.29.0, spglib 2.8.0, mace-torch 0.3.16, torch 2.14.1 CPU). The MACE model is `mace_mp(model="medium-mpa-0", default_dtype="float64")`. The MACE molecular-dynamics driver (`sse_md.py`) is not included; its trajectories are in NOMAD v2 under `mace_md/`.

2. DFT settings. We used Quantum ESPRESSO 7.5 with PBE-D3 and PAW pseudopotentials from pslibrary 1.0.0. Download the pseudopotentials from pslibrary and compare them against the file names and MD5 sums in `overnight/pseudo_md5.txt`. Each job's `.in` file states its cutoff and k-mesh.

3. Gate verdicts. Download the NOMAD v2 project files (FILES tab, download button) and unzip them into a folder `halide_models_v2/` next to this folder. Each row of its `manifest.csv` gives the dataset `path` of one job and the `source_job` path it had in this repository. Copy each `.in`/`.out` pair back to `source_job`:

   ```
   python -I - <<'EOF'
   import csv, shutil
   from pathlib import Path
   N = Path("halide_models_v2")  # adjust if unzipped elsewhere
   for r in csv.DictReader(open(N / "manifest.csv")):
       for ext in (".in", ".out"):
           src, dst = N / (r["path"] + ext), Path(r["source_job"] + ext)
           if src.exists():
               dst.parent.mkdir(parents=True, exist_ok=True)
               shutil.copy2(src, dst)
   EOF
   python -I discovery/scripts/verdict.py
   ```

   The restore skips `electrochemical_window/` and `reproduction/`, because `manifest.csv` does not list them and `verdict.py` does not read them. `verdict.py` rewrites `discovery/results/verdict.md` and `verdict.json` and ends with `self-check OK`. Gates G1 to G5 are defined in `README.md`.

4. Stability window. `discovery/scripts/window.py` reads the Na, Li and Br2 reference relaxations and the two gap SCFs from `overnight/m5_window/jobs/`. Copy them there from NOMAD v2, then run the script:
   `mkdir -p overnight/m5_window/jobs && cp halide_models_v2/electrochemical_window/*.in halide_models_v2/electrochemical_window/*.out overnight/m5_window/jobs/ && python -I discovery/scripts/window.py`
   The files in `reproduction/` are the same control SCF run on two laptops; no script reads them.

5. Paper numbers. After step 3 (and step 4 for the window values), run
   `python -I paper/claims.py && python -I paper/fill_preprint.py --draft`
   `claims.py` writes `paper/claims.json` and `claims.md`. The phonon fields stay empty unless the NOMAD v2 phonon summaries are copied in first: `dft_phonons/NaSbBr4/phonon_summary.json` to `anish_results/ionic_packet/anish_packet_run/S4/TP-NaSbBr4-01/`, and `dft_phonons/LiSbBr4/phonon_summary.json` to `.../S5/TP-LiSbBr4-01/`. `fill_preprint.py --draft` writes `paper/preprint_v2.draft.md` and marks any missing value as `[[key]]`. `python -I paper/figures/make_figures.py` redraws the figures.

6. New DFT jobs. The `overnight/make_*.py` scripts write inputs, and `overnight/run_night.py` runs a batch. Both need a QE build (`overnight/setup.sh`). `overnight/make_inputs.py` takes some starting cells from ML relaxations that are not included; the NOMAD v2 inputs carry the same geometries.

7. ML screening. Run `discovery/scripts/transplant_scan.py` and `screen_mace.py`. Scripts that query Materials Project need an API key in `~/.mp_api_key` (not included).

8. Cluster workflow (phonons, convergence checks). `anish_packet/go.sh` runs the steps listed in `anish_packet/PLAN` as Slurm jobs. Set the scratch and environment paths in `anish_packet/env.sh` first.

Status wording: all structures are predicted or proposed models (top two: TP-NaSbBr4-01 and TP-LiSbBr4-01). None has been synthesized or measured. No electrolyte is claimed.
