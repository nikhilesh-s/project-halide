"""Quick check that an unfinished relaxation resumes without the QE MPI_ABORT (2026-10-18 cell-noise fix).
  python3 -I overnight/test_resume.py m5_left/jobs/04_Na2SbBr6_agm005585966_ternary_vcrelax
Builds the resumed input exactly as run_night.py does (last clean geometry, cell noise removed), runs pw.x on 2 cores for
at most 3 minutes in a scratch folder, and prints PASS if the SCF starts (or FAIL with the last lines). Touches no job files."""
import os, re, subprocess, sys, tempfile, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_night import last_geometry, with_geometry

job = HERE / sys.argv[1]
inp, out = job.with_suffix(".in"), job.with_suffix(".out")
stale = lambda t: "The initial density is read from file" in t
cands = ([out] if out.exists() else []) + sorted(out.parent.glob(out.name + ".part*"), key=lambda f: int(f.name.rsplit(".part", 1)[1]), reverse=True)
geom = None
for f in cands:
    t = f.read_text(errors="replace")
    if not stale(t) and (geom := last_geometry(t)): print("geometry from", f.name); break
if not geom: sys.exit("no resumable geometry found")
tmp = Path(tempfile.mkdtemp(prefix="resume_test_"))
text = with_geometry(inp.read_text(), geom)
text = re.sub(r"prefix='[^']*'", "prefix='resume_test'", text).replace("__PSEUDO__", str(HERE / "pseudo")).replace("__TMP__", str(tmp)).replace("__MAXSEC__", "600")
print("cell used:\n" + re.search(r"CELL_PARAMETERS angstrom\n((?:.*\n){3})", text).group(1))
(tmp / "test.in").write_text(text)
log = tmp / "test.out"
p = subprocess.Popen([str(HERE / ".qe/bin/mpirun"), "-np", "2", str(HERE / ".qe/bin/pw.x"), "-i", str(tmp / "test.in")],
                     stdout=log.open("w"), stderr=subprocess.STDOUT, cwd=tmp, start_new_session=True)
t0 = time.time()
while time.time() - t0 < 180 and p.poll() is None:
    time.sleep(5)
    if "iteration #  1" in log.read_text(errors="replace"): break
ok = "iteration #  1" in log.read_text(errors="replace")
if p.poll() is None: os.killpg(p.pid, 9)  # stop mpirun AND its pw.x children
print("PASS: resumed geometry runs (got past the point where it used to abort and started the SCF); the queue will continue this job normally." if ok else
      "FAIL: still aborts. Last lines:\n" + "\n".join(log.read_text(errors="replace").splitlines()[-25:]))
