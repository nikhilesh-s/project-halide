"""One-command DFT verdict. Run from repo root:  discovery/.venv/bin/python -I discovery/scripts/verdict.py

Scans overnight/*/jobs/*.in|*.out (QE pw.x, PBE-D3), takes the last '!' energy (fresh-basis final SCF of a vc-relax),
and per chemical system prints every phase and, for every proposed model (job name contains "TP"):
  dE        vs the best listed (non-TP) structure of the same composition
  binaries  decomposition vs the computed binaries only
  competing decomposition vs the lowest balanced combination of ALL other-composition phases (linear programme)
  hull      vs everything else computed in the system = max(dE, competing)
All in meV/atom, negative = the proposed model is lower. A comparison only ever uses phases at the SAME setting
(ecutwfc + k-spacing tag); a reference that is absent at that setting prints "MISSING REFERENCE: <phase> at <cutoff> Ry".
Writes discovery/results/verdict.md and verdict.json. TP-* structures are PROPOSED MODELS; nothing here is a stability claim.
"""
import glob, importlib.util, json, re, sys
from collections import defaultdict
from pathlib import Path
from pymatgen.core import Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
RY = 13605.693123  # meV
HAL = ("F", "Cl", "Br", "I")
SKIP = {"test_air", "m5_window", "m5_phonons"}  # m5_phonons: 96-atom LiSbBr4 phonon displacements (not structures; would count as LiSbBr4 polymorphs).  test_air: build check on the Air, duplicate of night1/00. m5_window: elemental refs + gap SCFs for the stability window (overnight/make_window.py), analysed separately
_s = importlib.util.spec_from_file_location("make_inputs", ROOT / "overnight/make_inputs.py")
mi = importlib.util.module_from_spec(_s); _s.loader.exec_module(mi)  # reuse final_structure (last CELL/POSITIONS block)


def parse(inp):
    """One job -> dict. Status: done | not converged | running/partial | queued."""
    inp = Path(inp); out = inp.with_suffix(".out"); itxt = inp.read_text(); stem = inp.stem
    s = mi.final_structure(inp); scf = "'scf'" in itxt
    m = re.search(r"([A-Za-z0-9]+)_TP([A-Za-z]*)(\d+)", stem); a = re.search(r"agm\d+", stem); k = re.search(r"_k(\d)(\d\d)", stem)
    r = dict(job=f"{inp.parent.parent.name}/{stem}", formula=s.composition.reduced_formula, tp=bool(m),
             source=f"TP-{m.group(1)}-{m.group(2)}{m.group(3)}" if m else (a.group(0) if a else "(no id in job name)"),
             ecut=float(re.search(r"ecutwfc=([\d.]+)", itxt).group(1)), kspacing=f"0.{k.group(2)}" if k else "0.22",
             calc="scf" if scf else "vc-relax", natoms=len(s), status="queued", e_ry=None, e_ev_atom=None, sg=None, vol_atom=None)
    if out.exists():
        t = out.read_text(errors="replace")
        e = re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+) Ry", t, re.M)
        ec = re.search(r"kinetic-energy cutoff\s+=\s+([\d.]+)", t)
        if ec: r["ecut"] = float(ec.group(1))
        conv = ("convergence has been achieved" in t) if scf else ("bfgs converged" in t)
        paused = "stopped by user request" in t  # clean stop for a level switch; resumes from saved state
        stale = "The initial density is read from file" in t  # stale saved-state restart (2026-10-17 bug): never trusted
        r["status"] = "done" if ("JOB DONE" in t and conv and e and not paused and not stale) else ("not converged" if "JOB DONE" in t and not paused else "running/partial")
        if r["status"] == "done":
            if not scf: s = mi.final_structure(out)
            r.update(e_ry=float(e[-1]), e_ev_atom=float(e[-1]) * RY / 1000 / len(s), natoms=len(s),
                     sg=SpacegroupAnalyzer(s, 0.05).get_space_group_symbol(), vol_atom=s.volume / len(s))
    r["setting"] = f"{r['ecut']:.0f} Ry k{r['kspacing']}"
    return r


def decompose(target, phases):
    """Lowest-energy balanced combination of `phases` ({formula: E per atom, Ry}) giving one formula unit of `target`.
    Returns (energy per target atom, {formula: formula units}) or None when no balanced combination exists."""
    tc = Composition(target); names = sorted(phases); els = sorted({e.symbol for e in tc} | {e.symbol for n in names for e in Composition(n)})
    if not names: return None
    A = [[Composition(n)[el] for n in names] for el in els]
    res = linprog([phases[n] * Composition(n).num_atoms for n in names], A_eq=A, b_eq=[tc[el] for el in els], bounds=(0, None), method="highs")
    if not res.success: return None
    return res.fun / tc.num_atoms, {n: round(float(x), 4) for n, x in zip(names, res.x) if x > 1e-9}


def binaries_of(formula):
    c = Composition(formula); x = [e.symbol for e in c if e.symbol in HAL]; cat = [e.symbol for e in c if e.symbol not in HAL]
    if len(x) != 1 or len(cat) != 2: return []
    return [Composition(cat[0] + x[0]).reduced_formula, Composition(cat[1] + x[0] + "3").reduced_formula]  # AX, MX3 for A = Li/Na first


def analyse(rows):
    done = [r for r in rows if r["status"] == "done"]
    best = defaultdict(dict)  # setting -> formula -> lowest E/atom (Ry) among listed (non-TP) structures
    for r in done:
        if not r["tp"]:
            e = r["e_ry"] / r["natoms"]; b = best[r["setting"]]
            if e < b.get(r["formula"], (1e9,))[0]: b[r["formula"]] = (e, r["source"], r["job"])
    cands = defaultdict(list)
    for r in rows:
        if not r["tp"]: continue
        f, st = r["formula"], r["setting"]; sysels = {e.symbol for e in Composition(f)}; cut = f"{r['ecut']:.0f} Ry"
        v = dict(job=r["job"], setting=st, status=r["status"], formula=f)
        insys = {g: x for g, x in best[st].items() if {e.symbol for e in Composition(g)} <= sysels}
        pending = sorted({f"{q['formula']} ({q['job']}, {q['status']})" for q in rows if q["status"] != "done" and not q["tp"] and q["setting"] == st
                          and {e.symbol for e in Composition(q["formula"])} <= sysels
                          and not any(d["source"] == q["source"] and d["formula"] == q["formula"] and d["setting"] == st for d in done)})
        v["pending_in_system"] = pending
        # other proposed models of the same composition at this setting: a lower one would supersede this candidate
        v["rival_models"] = sorted(f"{q['source']} " + (f"{(q['e_ry'] / q['natoms'] - r['e_ry'] / r['natoms']) * RY:+.2f} vs this" if q["status"] == "done" and r["status"] == "done" else f"({q['status']}, {q['job']})")
                                   for q in rows if q["tp"] and q["formula"] == f and q["setting"] == st and q["source"] != r["source"])
        if r["status"] != "done":
            v.update(dE=None, binaries=None, competing=None, hull=None, note=f"candidate {r['status'].upper()} at {cut}")
        else:
            e = r["e_ry"] / r["natoms"]
            v["dE"] = round((e - insys[f][0]) * RY, 6) if f in insys else f"MISSING REFERENCE: listed {f} at {cut}"
            v["best_listed"] = insys[f][1] if f in insys else None
            others = {g: x[0] for g, x in insys.items() if g != f}
            bins = {g: x for g, x in others.items() if len(Composition(g)) == 2}
            for key, ph in (("binaries", bins), ("competing", others)):
                d = decompose(f, ph)
                if d: v[key] = round((e - d[0]) * RY, 6); v[key + "_route"] = d[1]
                else:
                    miss = [b for b in binaries_of(f) if b not in ph] or [f"balanced reference set for {f}"]
                    v[key] = "; ".join(f"MISSING REFERENCE: {b} at {cut}" for b in miss)
            v["ternaries_used"] = sorted(g for g in others if len(Composition(g)) == 3)
            num = [x for x in (v["dE"], v["competing"]) if isinstance(x, float)]
            v["hull"] = max(num) if len(num) == 2 else "MISSING REFERENCE (see dE / competing)"
        cands[r["source"]].append(v)
    return best, cands


def gates(vs):
    """4-gate line (gate definitions are in the paper, Section 2). G1 is the ML rattle gate (not recomputed: every TP job in a
    DFT queue passed it, gate_night3/4.md). G4 here = numerical convergence only; the literature scout is a manual check."""
    vs = sorted(vs, key=lambda v: (float(v["setting"].split()[0]), -float(v["setting"].split("k")[1]))); b = vs[0]; isn = lambda x: isinstance(x, float)  # base setting = lowest cutoff, coarsest k
    if b["status"] != "done": g2 = g3 = f"pending ({b['status']})"
    else:
        g2 = (f"pass ({b['dE']:+.2f} vs {b['best_listed']})" if b["dE"] < 0 else f"FAIL ({b['dE']:+.2f} vs {b['best_listed']})") if isn(b["dE"]) else b["dE"]
        if not isn(b["competing"]): g3 = b["competing"]
        elif b["competing"] >= 0: g3 = f"FAIL ({b['competing']:+.2f} vs {' + '.join(b['competing_route'])})"
        elif b["pending_in_system"]: g3 = f"pending ({b['competing']:+.2f} so far; not finished: {'; '.join(b['pending_in_system'])})"
        else: g3 = f"pass ({b['binaries']:+.2f} binaries; {b['competing']:+.2f} all competing; " + (f"ternaries tested: {', '.join(b['ternaries_used'])})" if b["ternaries_used"] else "NO competing ternary computed)")
    hi = [v for v in vs[1:] if v["status"] == "done" and isn(v["dE"]) and isn(v["competing"])]
    if hi and isn(b.get("dE")) and isn(b.get("competing")):
        h = hi[-1]; ok = h["dE"] < 0 and h["competing"] < 0 and abs(h["dE"] - b["dE"]) < 1 and abs(h["competing"] - b["competing"]) < 1
        g4 = f"{'converged pass' if ok else 'FAIL'} ({h['setting']}: dE {h['dE']:+.2f}, competing {h['competing']:+.2f}); scout: manual"
    else: g4 = f"pending (no second, tighter setting computed; base {b['setting']})"
    return dict(G1="pass (ML, not recomputed)", G2=g2, G3=g3, G4=g4, base_setting=b["setting"], rivals="; ".join(b["rival_models"]) or "none")


def main(write=True):
    # parked jobs that were started but not finished still count as unfinished competitors (2026-10-18: parking must not turn "pending" into "pass")
    started_parked = [f for f in glob.glob(str(ROOT / "overnight/parked/*/*.in")) if Path(f).with_suffix(".out").exists()]
    # jobs handed to helper machines (overnight/OUTSOURCED.txt) stay unfinished until their result is imported
    oz = ROOT / "overnight/OUTSOURCED.txt"
    started_parked += [str(ROOT / "overnight" / l.split("#")[0].strip()) for l in (oz.read_text().splitlines() if oz.exists() else [])
                       if l.strip() and not l.startswith("#")]
    started_parked = sorted(set(started_parked))
    rows = [parse(f) for f in sorted(glob.glob(str(ROOT / "overnight/*/jobs/*.in")) + started_parked) if Path(f).parent.parent.name not in SKIP
            and "vdw_corr='dft-d3'" in (t := Path(f).read_text()) and "dftd3_version" not in t]  # D3-BJ / no-D3 checks are another method: compared separately
    best, cands = analyse(rows)
    systems = sorted({"-".join(sorted(e.symbol for e in Composition(r["formula"]))) for r in rows if len(Composition(r["formula"])) == 3})
    fmt = lambda x: f"{x:+.2f}" if isinstance(x, float) else str(x)
    L = ["# DFT verdict (generated by `discovery/scripts/verdict.py`; do not edit by hand)", "",
         "QE PBE-D3, last `!` energy of each job. TP-* = PROPOSED MODEL. meV/atom, negative = proposed model lower. One setting per comparison.", ""]
    out = dict(systems={}, candidates={})
    for sy in systems:
        els = set(sy.split("-")); ph = [r for r in rows if {e.symbol for e in Composition(r["formula"])} <= els]
        L += [f"## {sy}", "", "| Phase | Source ID | Job | Setting | Atoms | Status | Space group | E (eV/atom) |", "|---|---|---|---|---|---|---|---|"]
        for r in sorted(ph, key=lambda r: (r["setting"], r["formula"], r["e_ev_atom"] or 0)):
            L.append(f"| {r['formula']} | {r['source']} | {r['job']} | {r['setting']} | {r['natoms']} | {r['status']} | {r['sg'] or ''} | {'' if r['e_ev_atom'] is None else format(r['e_ev_atom'], '.6f')} |")
        L.append(""); out["systems"][sy] = ph
        for c in sorted({r["source"] for r in ph if r["tp"]}):
            L += [f"**{c}**", "", "| Setting | Status | dE vs best listed | vs binaries | vs all competing phases | vs hull of everything else |", "|---|---|---|---|---|---|"]
            for v in sorted(cands[c], key=lambda v: v["setting"]):
                L.append(f"| {v['setting']} | {v['status']} | {fmt(v['dE'])}{' (' + v['best_listed'] + ')' if v.get('best_listed') else ''} | {fmt(v['binaries'])} | "
                         f"{fmt(v['competing'])}{' (' + ' + '.join(v['competing_route']) + ')' if v.get('competing_route') else ''} | {fmt(v['hull'])} |")
            g = gates(cands[c]); out["candidates"][c] = dict(rows=cands[c], gates=g)
            L += ["", f"Gates {c}: G1 {g['G1']} | G2 {g['G2']} | G3 {g['G3']} | G4 {g['G4']} | rival proposed models: {g['rivals']}", ""]
    L += ["## Gate summary", "", "| Candidate | Base setting | G1 | G2 | G3 | G4 | Rival proposed models (same composition) |", "|---|---|---|---|---|---|---|"]
    L += [f"| {c} | {d['gates']['base_setting']} | {d['gates']['G1']} | {d['gates']['G2']} | {d['gates']['G3']} | {d['gates']['G4']} | {d['gates']['rivals']} |" for c, d in sorted(out["candidates"].items())]
    if write:
        (ROOT / "discovery/results/verdict.md").write_text("\n".join(L) + "\n")
        (ROOT / "discovery/results/verdict.json").write_text(json.dumps(out, indent=1))
    return out, "\n".join(L)


def selfcheck(out):
    """E28 numbers: TP-NaBiCl4-01 at 60 Ry must give -19.82 (NaCl + BiCl3), -5.68 (agm088474086); competing was -12.79 (E28) and became -11.27 on 2026-10-19 when the P2_1/c Na3BiCl6 cell
    (agm003569174, m5_left/11) converged lower than the R3 cell: a real new competitor, not a regression."""
    assert decompose("NaBiCl4", {"NaCl": -1.0, "BiCl3": -1.0}) == (-1.0, {"BiCl3": 1.0, "NaCl": 1.0})   # LP balances 1:1
    assert decompose("NaBiCl4", {"NaCl": -1.0}) is None                                                # unbalanced -> missing
    v = next(v for v in out["candidates"]["TP-NaBiCl4-01"]["rows"] if v["setting"] == "60 Ry k0.22")
    for got, want in ((v["binaries"], -19.82), (v["competing"], -11.27), (v["dE"], -5.68)):
        assert abs(got - want) < 0.05, (got, want)
    assert v["best_listed"] == "agm088474086" and set(v["competing_route"]) == {"Na3BiCl6", "BiCl3"}, v
    print("self-check OK: TP-NaBiCl4-01 60 Ry", v["binaries"], v["competing"], v["dE"])


if __name__ == "__main__":
    o, md = main()
    print(md[md.index("## Gate summary"):])
    selfcheck(o)
    sys.exit(0)
