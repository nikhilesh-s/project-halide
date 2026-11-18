"""Electrochemical stability window + PBE gap of TP-NaSbBr4-01 / TP-LiSbBr4-01. Run from repo root:
  discovery/.venv/bin/python -I discovery/scripts/window.py [--selfcheck]
Grand-potential method (Zhu/He/Mo 2015; Richards/Mo 2016): open to M = Na/Li, mu_M = E(M metal) - V (V vs M/M+).
At each V, phi = E - n_M*mu_M; the candidate is stable if its phi equals the lowest balanced combination (LP over the
other elements) of all phases. Window = V range where that holds; products reported just outside it. If the candidate is
above hull at this setting, decomposition energy vs V is reported instead. All phases at ONE setting (80 Ry k0.22, PBE-D3):
verdict.py system rows + elemental refs from overnight/m5_window (01 Na, 02 Li, 03 Br2). Gap: 04/05 fixed-occupation SCF.
TP-* = PROPOSED MODEL; other TP models are not used as competitors. Prints 'PENDING: ...' for unfinished jobs, exit 0."""
import importlib.util, json, re, sys
from pathlib import Path
import numpy as np
from pymatgen.core import Composition
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
SET = "80 Ry k0.22"; TOL = 1e-5; VMAX = 6.0; STEP = 0.01  # TOL eV/f.u.; scan 0..VMAX V
J = ROOT / "overnight/m5_window/jobs"
CAND = {"Na": ("TP-NaSbBr4-01", "Br-Na-Sb", "01_Na_bcc_ref_vcrelax", "04_NaSbBr4_TP01_gap_scf"),
        "Li": ("TP-LiSbBr4-01", "Br-Li-Sb", "02_Li_bcc_ref_vcrelax", "05_LiSbBr4_TP01_gap_scf")}
BR2 = "03_Br2_Cmca_ref_vcrelax"


def lp(phases, target, open_el, mu):
    """phases [(label, Composition, E_fu eV)] -> (min grand potential per target f.u., {label: f.u.})."""
    closed = sorted({e.symbol for e in target if e.symbol != open_el})
    A = [[c[el] for _, c, _ in phases] for el in closed]
    cost = [e - c[open_el] * mu for _, c, e in phases]
    r = linprog(cost, A_eq=A, b_eq=[target[el] for el in closed], bounds=(0, None), method="highs")
    assert r.success, "no balanced decomposition"
    return r.fun, {l: round(float(x), 4) for (l, _, _), x in zip(phases, r.x) if x > 1e-9}


def window(phases, ci, open_el, e_open):
    """phases: list incl. candidate at index ci; e_open = E per atom of open-element metal. Returns dict."""
    _, tc, te = phases[ci]
    g = lambda V: lp(phases, tc, open_el, e_open - V)[0] - (te - tc[open_el] * (e_open - V))  # <= 0; 0 = stable
    stab = lambda V: g(V) >= -TOL
    Vs = np.arange(0, VMAX + 1e-9, STEP); s = [stab(V) for V in Vs]
    out = dict(candidate=phases[ci][0], open=open_el)
    if not any(s):
        out.update(stable=False, max_stability_V=float(Vs[int(np.argmax([g(V) for V in Vs]))]),
                   decomp_vs_V={f"{V:.1f}": round(-g(V) / tc.num_atoms * 1000, 2) for V in np.arange(0, VMAX + 1e-9, 0.5)},
                   note="above hull at this setting: meV/atom above grand-potential hull vs V (V vs M/M+)")
        return out
    i0 = s.index(True); i1 = len(s) - 1 - s[::-1].index(True)
    def bis(a, b):  # a unstable, b stable (or vice versa); find edge to 1e-4 V
        sa = stab(a)
        for _ in range(30):
            m = (a + b) / 2
            if stab(m) == sa: a = m
            else: b = m
        return (a + b) / 2
    lo = 0.0 if i0 == 0 else bis(Vs[i0 - 1], Vs[i0]); hi = VMAX if i1 == len(Vs) - 1 else bis(Vs[i1], Vs[i1 + 1])
    d = 0.01; prod = lambda V: lp(phases, tc, open_el, e_open - V)[1]
    out.update(stable=True, V_reduction=round(lo, 3), V_oxidation=round(hi, 3), window_V=round(hi - lo, 3),
               reduction_products=None if i0 == 0 else prod(lo - d), oxidation_products=None if i1 == len(Vs) - 1 else prod(hi + d))
    if i0 == 0: out["note_red"] = "stable down to 0 V (vs metal): no reduction limit above the metal"
    if i1 == len(Vs) - 1: out["note_ox"] = f"stable to the scan limit {VMAX} V"
    return out


def gap(name):
    p = J / f"{name}.out"
    if not p.exists(): return None
    t = p.read_text(errors="replace"); m = re.findall(r"highest occupied, lowest unoccupied level \(ev\):\s+(-?[\d.]+)\s+(-?[\d.]+)", t)
    if "JOB DONE" not in t or "convergence has been achieved" not in t or not m: return None
    return dict(homo=float(m[-1][0]), lumo=float(m[-1][1]), gap_eV=round(float(m[-1][1]) - float(m[-1][0]), 3))


def selfcheck():
    C = Composition  # A = open element, B closed. Analytic: AB2+3A->2A2B at mu=-1/3 (V=1/3); AB2->1/2 AB4+1/2 A at mu=-1 (V=1)
    ph = [("AB2", C("AB2"), -3.0), ("A", C("A"), 0.0), ("B", C("B"), 0.0), ("A2B", C("A2B"), -2.0), ("AB4", C("AB4"), -5.0)]
    w = window(ph, 0, "A", 0.0)
    assert w["stable"] and abs(w["V_reduction"] - 1 / 3) < 1e-3 and abs(w["V_oxidation"] - 1.0) < 1e-3, w
    assert "A2B" in w["reduction_products"] and "AB4" in w["oxidation_products"], w
    w2 = window([("AB2", C("AB2"), -1.0)] + ph[1:], 0, "A", 0.0); assert not w2["stable"]  # raised above hull -> reports decomposition instead
    print("selfcheck OK:", w["V_reduction"], w["V_oxidation"])


def main():
    spec = importlib.util.spec_from_file_location("verdict", ROOT / "discovery/scripts/verdict.py"); V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
    refs = {n: V.parse(J / f"{n}.in") for n in (CAND["Na"][2], CAND["Li"][2], BR2)}
    pending = [f"{n} (vc-relax {r['status']})" for n, r in refs.items() if r["status"] != "done"] + [f"{c[3]} (SCF not done)" for c in CAND.values() if gap(c[3]) is None]
    systems, _ = V.main(write=False); systems = systems["systems"]; res = {}; lines = ["# Stability window and PBE gap (generated by `discovery/scripts/window.py`)", "",
        f"PBE-D3, setting {SET}. TP-* = PROPOSED MODEL; voltages vs M/M+ from a 0 K grand-potential hull (no kinetics, no interphase passivation).", ""]
    for M, (cid, sy, mref, gjob) in CAND.items():
        r = dict(setting=SET)
        if refs[mref]["status"] == "done" and refs[BR2]["status"] == "done":
            ph = []
            for row in systems[sy]:
                if row["status"] != "done" or row["setting"] != SET or (row["tp"] and row["source"] != cid): continue
                c = Composition(row["formula"]); ph.append((f"{row['formula']} [{row['source']}]", c, row["e_ry"] / row["natoms"] * V.RY / 1000 * c.num_atoms))
            for n in (mref, BR2):
                x = refs[n]; assert x["setting"] == SET or x["ecut"] == 80.0, x["setting"]; c = Composition(x["formula"]); ph.append((f"{x['formula']} [{n[:2]}]", c, x["e_ry"] / x["natoms"] * V.RY / 1000 * c.num_atoms))
            keep = {}  # lowest E per (formula, TP flag): duplicates of the same structure type (e.g. parked copies) collapse
            for l, c, e in ph:
                k = l if "TP-" in l else c.reduced_formula + ("*" if c.reduced_formula in ("Sb", "Br", M) else "")
                if k not in keep or e / c.num_atoms < keep[k][2] / keep[k][1].num_atoms: keep[k] = (l, c, e)
            ph = list(keep.values()); r["phases"] = [p[0] for p in ph]
            ci = next(i for i, p in enumerate(ph) if cid in p[0]); emet = next(e / c.num_atoms for l, c, e in ph if l.startswith(M + " ["))
            r["window"] = window(ph, ci, M, emet)
        g = gap(gjob)
        if g: r["pbe_gap"] = g
        if r.get("window") or g: res[cid] = r
    for cid, r in res.items():
        lines += [f"## {cid}", ""]
        if "window" in r:
            w = r["window"]; lines.append(f"Phases used ({len(r['phases'])}): {', '.join(r['phases'])}")
            lines.append(f"- Window vs {w['open']}/{w['open']}+: {w['V_reduction']} to {w['V_oxidation']} V ({w['window_V']} V wide)" if w["stable"] else f"- ABOVE HULL at this setting; {w['note']}: {w['decomp_vs_V']}")
            for k in ("reduction_products", "oxidation_products", "note_red", "note_ox"):
                if w.get(k): lines.append(f"- {k.replace('_', ' ')}: {w[k]}")
        if "pbe_gap" in r: lines.append(f"- PBE gap (fixed occupation, D3): {r['pbe_gap']['gap_eV']} eV (underestimates the true gap)")
        lines.append("")
    if pending: lines += ["PENDING: " + "; ".join(pending)]; print("PENDING:", "; ".join(pending))
    if res:
        (ROOT / "discovery/results/window.json").write_text(json.dumps(dict(results=res, pending=pending), indent=1)); (ROOT / "discovery/results/window.md").write_text("\n".join(lines) + "\n"); print("\n".join(lines))


if __name__ == "__main__":
    selfcheck() if "--selfcheck" in sys.argv else main()
    sys.exit(0)
