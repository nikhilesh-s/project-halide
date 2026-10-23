#!/usr/bin/env python3
"""Run queued pw.x jobs until a wall-clock stop time, then exit. Stdlib only, works offline.

  python3 -I overnight/run_night.py night1 --stop 07:15

Each job is skipped if its .out already says JOB DONE. An unfinished vc-relax is resumed from the last
geometry it printed. pw.x is told to stop itself (max_seconds) before the deadline; a hard kill is the backstop.
"""
import argparse, datetime as dt, json, os, re, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
RY_TO_MEV = 13605.693123
ANISH_NaBiCl4_D3_RY = -2667.12492286  # last "!" (final SCF on the relaxed cell, fresh basis) in dft_nabicl4/NaBiCl4_d3.out, 12 atoms.
# Note: his "Final enthalpy" is -2667.13299844 Ry; the 9 meV/atom gap is the basis-set change after vc-relax at 60 Ry.

def seconds_until(hhmm):
    h, m = map(int, hhmm.split(":"))
    now = dt.datetime.now()
    stop = now.replace(hour=h, minute=m, second=0, microsecond=0)
    if stop <= now:
        stop += dt.timedelta(days=1)
    return (stop - now).total_seconds()

def last_geometry(out_text):
    """Return (cell_block, positions_block) text from the last relaxation step, or None."""
    if "CELL_PARAMETERS" not in out_text:
        return None
    blk = out_text[out_text.rindex("CELL_PARAMETERS"):].splitlines()
    m = re.search(r"alat=\s*([\d.]+)", blk[0])
    scale = float(m.group(1)) * 0.529177210903 if m else 1.0
    cell = ["  " + " ".join(f"{float(x) * scale:.9f}" for x in l.split()) for l in blk[1:4]]
    try:
        i = next(n for n, l in enumerate(blk) if l.startswith("ATOMIC_POSITIONS"))
    except StopIteration:
        return None
    pos = []
    for l in blk[i + 1:]:
        w = l.split()
        if len(w) != 4 or not w[0].isalpha():
            break
        pos.append("  " + l.strip())
    return "\n".join(cell), "\n".join(pos)

def clean_cell(cell):
    """Zero cell components below 1e-5 A. Resumed cells carry BFGS noise like 0.000000989 in off-diagonals; QE then finds a
    slightly different symmetry on different MPI ranks and aborts right after the RAM estimate (m5_left/04 Na2SbBr6, 2026-10-18;
    night2 10a1 earlier). Real cell components are > 0.1 A, so this only removes the noise."""
    out = []
    for line in cell.splitlines():
        v = line.split()
        out.append("  " + " ".join(f"{0.0 if abs(float(x)) < 1e-5 else float(x):.10f}" for x in v) if len(v) == 3 else line)
    return "\n".join(out)

def with_geometry(inp, geom):
    cell, pos = geom
    cell = clean_cell(cell)
    inp = re.sub(r"(CELL_PARAMETERS angstrom\n)(?:.*\n){3}", lambda m: m.group(1) + cell + "\n", inp)
    return re.sub(r"(ATOMIC_POSITIONS crystal\n)(?:  \w+ .*\n)+", lambda m: m.group(1) + pos + "\n", inp)

def final_energy_ry(out_text):
    e = re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", out_text, re.M)
    return float(e[-1]) if e else None

def summarize(jobdir):
    rows = {}
    for out in sorted(jobdir.glob("*.out")):
        t = out.read_text(errors="replace")
        nat = re.search(r"number of atoms/cell\s+=\s+(\d+)", t)
        e = final_energy_ry(t)
        rows[out.stem] = {
            "job_done": "JOB DONE" in t,
            "relax_converged": "bfgs converged" in t or "calculation='scf'" in (jobdir / (out.stem + ".in")).read_text(),
            "natoms": int(nat.group(1)) if nat else None,
            "final_energy_Ry": e,
            "energy_per_atom_Ry": None if not (e and nat) else e / int(nat.group(1)),
            "energy_per_atom_meV_vs_Anish_NaBiCl4_d3": None if not (e and nat and "_NaBiCl4_" in out.stem) else round(
                (e / int(nat.group(1)) - ANISH_NaBiCl4_D3_RY / 12) * RY_TO_MEV, 3),
            "wall": (re.findall(r"PWSCF\s+:.*?CPU\s+(.*?) WALL", t) or [None])[-1],
        }
    (jobdir.parent / "summary.json").write_text(json.dumps(rows, indent=2))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("night")
    ap.add_argument("--stop", default="07:15", help="local wall-clock time HH:MM to be finished by")
    ap.add_argument("--np", type=int, default=0, help="MPI processes (default: performance cores)")
    a = ap.parse_args()
    jobdir = HERE / a.night / "jobs"
    pw = HERE / ".qe/bin/pw.x"
    mpirun = HERE / ".qe/bin/mpirun"
    assert pw.exists() and mpirun.exists(), "run overnight/setup.sh first"
    n = a.np or int(subprocess.run(["sysctl", "-n", "hw.perflevel0.physicalcpu"], capture_output=True, text=True).stdout or 4)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PATH=f"{HERE / '.qe/bin'}:{os.environ['PATH']}")
    tmp = HERE / "tmp"; tmp.mkdir(exist_ok=True)
    deadline = time.time() + seconds_until(a.stop)
    for inp in sorted(jobdir.glob("*.in")):
        out = inp.with_suffix(".out")
        prev = out.read_text(errors="replace") if out.exists() else ""
        if "JOB DONE" in prev and "stopped by user request" not in prev and "The initial density is read from file" not in prev and (
                "bfgs converged" in prev or ("'scf'" in inp.read_text() and "convergence has been achieved" in prev)):
            print(f"skip (done): {inp.name}"); continue
        left = deadline - time.time()
        if left < 900:
            print(f"stop: only {left:.0f}s left before {a.stop}"); break
        np_job = n
        if "phonon" in inp.name:  # 96-atom LiSbBr4 single-point, cannot resume. RAM ~ 3.2 GB/proc + 5.6 GB (a volunteer 12 procs 44 GB, a volunteer 10 procs 38 GB)
            mem_gb = int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout or 0) / 2**30
            np_job = min(n, int((mem_gb - 5 - 5.6) / 3.2))  # keep ~5 GB for macOS; 24 GB M5 -> 4 procs (~18 GB)
            need_h = 3.9 * 12 / max(np_job, 1) * 1.3  # 3h54m on a volunteer's 12 procs, scaled, 30% margin
            print(f"phonon job: RAM {mem_gb:.0f} GB -> np {np_job}, needs ~{need_h:.0f} h", flush=True)
            if np_job < 3:
                print(f"skip: {inp.name} needs >= ~20 GB RAM"); continue
            if left < need_h * 3600:
                print(f"skip for now: {inp.name} needs ~{need_h:.0f} h, only {left/3600:.1f} h left"); continue
        text = inp.read_text()
        flag = tmp / (inp.stem + ".graceful")  # written by the dashboard after a clean stop (prefix.EXIT): QE saved its state
        # 2026-10-17: restart_mode='restart' re-read a saved charge density for a different cell (air/22 Na3BiCl6 converged
        # 2.5 eV/atom too high). Never restart from saved state; resume from the last geometry of a segment that started clean.
        stale = lambda t: "The initial density is read from file" in t
        geom = last_geometry(prev) if prev and not stale(prev) else None
        if not geom:  # a failed or stale restart: fall back to the newest earlier clean part
            for part in sorted(out.parent.glob(out.name + ".part*"), key=lambda f: int(f.name.rsplit(".part", 1)[1]), reverse=True):  # suffix = time written; mtime resets on git pull
                t = part.read_text(errors="replace")
                geom = None if stale(t) else last_geometry(t)
                if geom: break
        flag.unlink(missing_ok=True)
        if geom and "vc-relax" in text:
            text = with_geometry(text, geom)
            out.rename(out.with_suffix(f".out.part{int(time.time())}"))
            print(f"resume from last geometry: {inp.name}")
        text = (text.replace("__PSEUDO__", str(HERE / "pseudo")).replace("__TMP__", str(tmp))
                    .replace("__MAXSEC__", str(int(left - 600))))
        run_in = tmp / inp.name
        run_in.write_text(text)
        print(f"{dt.datetime.now():%H:%M} start {inp.name} np={np_job} budget={left / 3600:.1f}h", flush=True)
        with open(out, "w") as fo:
            try:
                subprocess.run(([str(mpirun), "-np", str(np_job), str(pw), "-nk", "2" if np_job % 2 == 0 else "1"] if np_job > 1 else [str(pw)]) + ["-i", str(run_in)],  # np 1: no mpirun, so no network needed
                               stdout=fo, stderr=subprocess.STDOUT, env=env, timeout=max(60, left - 120))
            except subprocess.TimeoutExpired:
                print("hard stop at deadline")
    rows = summarize(jobdir)
    print(json.dumps(rows, indent=2))

if __name__ == "__main__":
    sys.exit(main())
