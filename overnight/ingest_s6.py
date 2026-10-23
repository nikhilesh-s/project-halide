"""File Anish's S6 results (disp-NNN.in/.out from the tarball of `go.sh pack`) under overnight/s6/jobs/ so verdict.py reads them.
Run from repo root:  discovery/.venv/bin/python -I overnight/ingest_s6.py <extracted_run_dir_or_S6_top2_dir> [--dry-run]
Rules (docs/review_oct18/bottom_up.md A.3 #1-3, #5):
  * disp-NNN -> laptop job by anish_packet/S6/top2/MANIFEST.tsv, then VERIFIED by content: the cluster .in must equal the packet
    input once pseudo_dir / outdir / prefix / max_seconds are set aside. Mismatch = ERROR, nothing filed.
  * done = verdict.py's rule (JOB DONE + bfgs/scf converged + a '!' energy + no stale-density restart). A max_seconds stop prints
    JOB DONE too: flagged "timed out", copied to s6/requeue/ (last geometry kept) for the Air, never filed as a result.
  * Done ones go to s6/jobs/<manifest stem>.{in,out}; .in = packet input with laptop placeholders. Writes ONLY under overnight/s6/.
  * Duplicates: the same structure under another name (number prefix stripped, e.g. new/SbBr3_X <-> air/21a_SbBr3_X) done on the
    Air/M5 wins: the cluster copy goes to s6/reproduction/ (not read by verdict.py) and the meV/atom difference is reported.
    Otherwise the cluster copy is filed in s6/jobs. Re-running moves files to wherever the rule says now, so it is idempotent.
  * 'listed' jobs: StructureMatcher of the relaxed structure vs the TP model's relaxed structure.
"""
import argparse, importlib.util, re, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent; O = ROOT / "overnight"; S6 = O / "s6"
RY = 13605.693123  # meV
TP = {"NaSbBr4": ("TP-NaSbBr4-01", "night3/jobs/09_NaSbBr4_TP01_hit_vcrelax.out"), "LiSbBr4": ("TP-LiSbBr4-01", "night3/jobs/02_LiSbBr4_TP01_hit_vcrelax.out")}
_s = importlib.util.spec_from_file_location("make_inputs", O / "make_inputs.py"); mi = importlib.util.module_from_spec(_s); _s.loader.exec_module(mi)

norm = lambda t: re.sub(r"(pseudo_dir|outdir|prefix|max_seconds)\s*=\s*\S+", r"\1=X", t).strip()
base = lambda stem: re.sub(r"^\d+[a-z]?_", "", stem)  # 21a_SbBr3_x_vcrelax == new/SbBr3_x_vcrelax


def status(t, scf):
    """Same rule as verdict.parse(); 'timed out' = QE hit max_seconds but still printed JOB DONE."""
    e = re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", t, re.M)
    conv = ("convergence has been achieved" in t) if scf else ("bfgs converged" in t)
    if "JOB DONE" not in t or "stopped by user request" in t: return "partial", None
    if conv and e and "The initial density is read from file" not in t: return "done", float(e[-1])
    return ("timed out" if "Maximum CPU time exceeded" in t else "not converged"), None


def find_run(p):
    p = Path(p)
    for c in [p, p / "S6/top2", *sorted(p.glob("*/S6/top2")), *sorted(p.glob("*/*/S6/top2"))]:
        if (c / "disp-000.out").exists() or list(c.glob("disp-*.out")): return c
    sys.exit(f"no disp-*.out under {p}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--dry-run", action="store_true"); a = ap.parse_args()
    run = find_run(a.run); pk = ROOT / "anish_packet/S6/top2"
    man = dict(l.split("\t") for l in (pk / "MANIFEST.tsv").read_text().split("\n")[1:] if l)
    others = {}  # base name -> [(job, in, out)] for every non-s6 job, including started parked ones
    for f in sorted(O.glob("*/jobs/*.in")) + sorted(O.glob("parked/*/*.in")):
        if f.parent.parent.name != "s6": others.setdefault(base(f.stem), []).append(f)
    rows, errs, requeue, move = [], [], [], []
    for disp, src in man.items():
        o, i = run / f"{disp}.out", run / f"{disp}.in"; stem = Path(src).name; r = dict(disp=disp, job=src, st="MISSING", note="")
        rows.append(r)
        if not o.exists(): continue
        pin = (pk / f"{disp}.in").read_text()
        if not i.exists() or norm(i.read_text()) != norm(pin):
            r.update(st="ERROR", note="cluster .in does not match the packet input" if i.exists() else "cluster .in missing"); errs.append(disp); continue
        t = o.read_text(errors="replace"); scf = "calculation='scf'" in pin
        nat = re.search(r"nat=(\d+)", pin).group(1); onat = re.search(r"number of atoms/cell\s+=\s+(\d+)", t)
        if not onat or onat.group(1) != nat:
            r.update(st="ERROR", note=f".out has {onat and onat.group(1)} atoms, input {nat}"); errs.append(disp); continue
        st, e = status(t, scf); r["st"] = st; r["nat"] = int(nat)
        if st == "timed out": requeue.append(src); move.append((o, S6 / "requeue" / f"{stem}.out"))
        if st != "done": continue
        r["e"] = e
        # the Air/M5 twin(s): same structure job under another name; a done twin wins
        twin = None
        for f in others.get(base(stem), []):
            ts, te = status(f.with_suffix(".out").read_text(errors="replace"), scf) if f.with_suffix(".out").exists() else ("queued", None)
            if ts == "done": twin = (f, te)
        if twin:
            d = (e - twin[1]) / int(nat) * RY
            r.update(st="done (twin done)", note=f"{twin[0].relative_to(O)}: {d:+.3f} meV/atom" + ("" if abs(d) < 0.5 else "  ** >0.5, look at it **"))
            dest = S6 / "reproduction"
        else: dest = S6 / "jobs"
        lap = re.sub(r"prefix='[^']*'", f"prefix='{stem}'", pin); lap = re.sub(r"outdir='[^']*'", f"outdir='__TMP__/{stem}'", lap)
        lap = re.sub(r"max_seconds=\d+", "max_seconds=__MAXSEC__", lap)
        move += [(lap, dest / f"{stem}.in"), (o, dest / f"{stem}.out")]
        if "listed" in stem and dest == S6 / "jobs" or "listed" in stem and twin is None:
            pass
        if "listed" in stem and not scf:
            f = next(f for f in TP if f in stem); name, ref = TP[f]
            from pymatgen.analysis.structure_matcher import StructureMatcher
            ok = StructureMatcher().fit(mi.final_structure(o), mi.final_structure(O / ref))
            r["note"] += f" {name}: {'MATCH' if ok else 'DISTINCT'}"
            print(("!" * 70 + f"\n!! {src} vs {name}: " + ("MATCH: relaxed into the TP structure, TP is a LISTED structure" if ok else "DISTINCT") + "\n" + "!" * 70))
    for s, d in move:  # write only under overnight/s6; a file filed earlier under the other folder is moved, not duplicated
        if a.dry_run: continue
        d.parent.mkdir(parents=True, exist_ok=True); stem = d.stem if d.suffix == ".in" else d.name
        for sib in ("jobs", "reproduction"):
            if S6 / sib != d.parent and (S6 / sib / d.name).exists() and d.parent.name != "requeue": (S6 / sib / d.name).unlink()
        d.write_text(s) if isinstance(s, str) else shutil.copyfile(s, d)
    L = ["| disp | laptop job | status | note |", "|---|---|---|---|"] + [f"| {r['disp']} | {r['job']} | {r['st']} | {r['note'].strip()} |" for r in rows]
    cnt = {k: sum(r["st"].startswith(k) for r in rows) for k in ("done", "timed out", "not converged", "partial", "MISSING", "ERROR")}
    head = f"{len(rows)} MANIFEST rows: " + ", ".join(f"{v} {k}" for k, v in cnt.items() if v)
    print(head + (" (dry run, nothing written)" if a.dry_run else "")); print("\n".join(L))
    if requeue: print("\nREQUEUE ON THE AIR (timed out, last geometry in overnight/s6/requeue/):\n  " + "\n  ".join(requeue))
    if not a.dry_run:
        S6.mkdir(exist_ok=True)
        (S6 / "INGEST_LOG.md").write_text(f"# S6 ingest log (generated by `overnight/ingest_s6.py`; source `{run.name}`)\n\n{head}\n\n" + "\n".join(L) + "\n" +
            ("\nRequeue on the Air (timed out): " + ", ".join(requeue) + "\n" if requeue else ""))
    sys.exit(1 if errs else 0)


if __name__ == "__main__":
    main()
