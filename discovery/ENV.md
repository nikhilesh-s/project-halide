# Discovery environment

- Python 3.11.14 (uv venv at discovery/.venv); run scripts with `discovery/.venv/bin/python -I`
- pymatgen 2026.9.24, ase 3.29.0, spglib 2.8.0, phonopy 4.8.3, mace-torch 0.3.16, torch 2.14.1 (CPU), numpy 2.4.6
- MACE `mace_mp(model="medium-mpa-0", device="cpu", default_dtype="float64")` loads OK (cached at ~/.cache/mace/macempa0mediummodel)
- Note: the `uv` on PATH may block `uv pip`; installed via /path/to/home/.local/bin/uv pip install.
