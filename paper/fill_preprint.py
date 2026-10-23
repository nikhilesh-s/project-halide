"""Fill {{key}} placeholders in paper/preprint_v2.md from paper/claims.json and paper/figures/data only.
Run from the repo root:
    discovery/.venv/bin/python -I paper/fill_preprint.py           -> paper/preprint_v2.filled.md, exit 1 on any missing/PENDING key
    discovery/.venv/bin/python -I paper/fill_preprint.py --draft   -> paper/preprint_v2.draft.md (missing keys shown as [[key]])
Placeholder syntax: {{key}} or {{key|fmt}}; fmt is a Python format spec, optionally prefixed by "abs" (e.g. {{na_comp|abs.1f}}).
Key kinds:
  alias        (table ALIAS below)            -> a field of claims.json "claims", or a value derived from figures/data
  T1.<row>.<col>                              -> paper/figures/data/T1_method_error.csv (rows T1_ROWS, cols zero|bj|none)
  meta.<NAME>, search.<name>, control.<a.b>  -> top-level blocks "meta" / "search" / "control" of claims.json (claims.py must emit them)
  struct.<a.b>                                -> paper/figures/data/T2_structures.json (written by the figures/structure step)
  include:<file>                              -> whole file paper/figures/data/<file> pasted in (first "# git" stamp line dropped)
  window_<F>, gap_<F>, win_lo/hi/width_<F>    -> discovery/results/window.json (F = NaSbBr4 | LiSbBr4); PENDING if absent
Conditional text (may nest): {{?key}}A{{:}}B{{/?}} renders A if key has a value (not 0/""/"none"), else B;
{{?key=value}}A{{:}}B{{/?}} renders A if the value equals `value`. An =-condition whose key is missing is reported as MISSING
(the choice must be explicit in claims.json); a bare condition just takes the B branch.
A key is PENDING if its value says pending, or it comes from a claim whose status is PENDING."""
import argparse, csv, json, re, sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

P = Path(__file__).resolve().parent
J = json.loads((P / "claims.json").read_text())
C = J["claims"]
NA, LI = "TP-NaSbBr4-01", "TP-LiSbBr4-01"


def table(name):
    lines = [l for l in (P / "figures/data" / name).read_text().splitlines() if not l.startswith("#")]
    return list(csv.DictReader(lines))


T1 = {r["quantity (meV/atom)"]: r for r in table("T1_method_error.csv")}
F1 = table("F1_competing_phases.csv")
T1_ROWS = {"na_bin": "TP-NaSbBr4-01 vs NaBr + lowest computed SbBr3",
           "na_listed": "TP-NaSbBr4-01 vs lowest listed NaSbBr4 (agm003733885)",
           "li_bin": "TP-LiSbBr4-01 vs LiBr + SbBr3 (lowest LiBr computed in each column)",
           "li_rs": "TP-LiSbBr4-01 vs LiBr rocksalt (agm003157564) + SbBr3",
           "li_wz": "TP-LiSbBr4-01 vs LiBr wurtzite (agm002091979) + SbBr3",
           "li_listed": "TP-LiSbBr4-01 vs lowest listed LiSbBr4 (agm005739369)",
           "sbbr3_order": "SbBr3: Pnma (expt, agm003218591) minus R-3 (agm002139075)",
           "libr_order": "LiBr: rocksalt (expt, agm003157564) minus wurtzite (agm002091979)"}
T1_COLS = {"zero": "PBE-D3 (zero damping)", "bj": "PBE-D3(BJ)", "none": "PBE (no D3)"}
PENDING, MISSING = "PENDING", "MISSING"
_T2 = P / "figures/data/T2_structures.json"
BLOCKS = {"meta": J.get("meta", {}), "search": J.get("search", {}), "control": J.get("control", {}), "nonlocal": J.get("nonlocal_check", {}),
          "struct": json.loads(_T2.read_text()) if _T2.exists() else {}}


def claim(tp, *path):
    v = C[tp]
    for k in path:
        if not isinstance(v, dict) or k not in v:
            return None, MISSING
        v = v[k]
    if v is None or (isinstance(v, str) and v.upper().startswith(("PENDING", "MISSING"))):
        return None, PENDING if v is None or v.upper().startswith("PENDING") else MISSING
    return v, PENDING if C[tp]["status"] != "FINAL" else None


def f1(system, source):
    r = next((r for r in F1 if r["system"] == system and r["source"] == source), None)
    if r is None: return None, MISSING
    if r["status"] != "done": return None, PENDING
    return float(r["dE_vs_binaries_meV_atom"]), None


def na_ternaries(fn):
    v = [float(r["dE_vs_binaries_meV_atom"]) for r in F1 if r["system"] == "Na-Sb-Br" and r["status"] == "done"
         and r["role"] != "TP" and r["formula"] not in ("NaBr", "SbBr3")]
    return (fn(v), None) if v else (None, MISSING)


def li_listed2_gap():  # distinct listed C2/c LiSbBr4 (agm003740876) minus the TP, same zero
    (a, ea), (b, eb) = f1("Li-Sb-Br", "agm003740876"), f1("Li-Sb-Br", "TP-LiSbBr4-01")
    return (a - b, None) if not (ea or eb) else (None, ea or eb)


def li_phonon_status():
    ph = C[LI]["phonons"]
    if not ph: return "not complete at the time of writing", PENDING
    if ph.get("verdict") != "PASS":  # never print "no instability" for a FAIL: stop the render so the text is rewritten by hand
        return f"imaginary modes found (mesh minimum {ph['min_mesh_THz']:+.3f} THz)", "PHONON VERDICT " + str(ph.get("verdict"))
    return (f"no instability found at {'x'.join(map(str, ph['supercell']))} sampling (mesh minimum "
            f"{ph['min_mesh_THz']:+.3f} THz)", None if C[LI]["status"] == "FINAL" else PENDING)


def li_ph(field):  # LiSbBr4 phonon numbers; same verdict guard as li_phonon_status (text says "no imaginary mode")
    _, err = li_phonon_status()
    if err: return None, err
    return claim(LI, "phonons", field)


def pending_ids(tp, formula=None):  # {agm id: formula} of distinct pending phases; one phase is often queued on two machines
    v, e = claim(tp, "pending")
    if v is None: return None, e
    ids = {}
    for p in v:
        f, m = p.split(" ")[0], re.search(r"agm\d+", p)
        if formula in (None, f): ids[m.group(0) if m else p] = f
    return ids, e


def sub_formula(f): return re.sub(r"(\d+)", r"~\1~", f)


def pending_names(tp):  # "3 SbBr~3~ (agm..., agm..., agm...); 2 NaBr (...)"
    ids, e = pending_ids(tp)
    if ids is None: return None, e
    by = {}
    for i, f in ids.items(): by.setdefault(f, []).append(i)
    return "; ".join(f"{len(v)} {sub_formula(f)} ({', '.join(sorted(v))})" for f, v in sorted(by.items())) or "none", e


def n_pending(formula=None):  # distinct pending phases over both systems
    ids, err = {}, None
    for tp in (NA, LI):
        v, e = pending_ids(tp, formula); err = err or e
        if v: ids |= v
    return len(ids), err


def sbbr3_flip(tp):  # how far a new SbBr3 polymorph must lie below the reference to flip the AX + SbBr3 margin (SbBr3 = 4 of 6 atoms)
    v, e = claim(tp, "vs_competing")
    return (abs(v) * 6 / 4, e) if v is not None else (v, e)


_WIN = P.parent / "discovery/results/window.json"  # written by discovery/scripts/window.py; absent -> PENDING
WIN = json.loads(_WIN.read_text()).get("results", {}) if _WIN.exists() else {}


def products(d):  # {"NaBr [agm...]": f.u., ...} -> "NaBr + Sb"
    return " + ".join(sub_formula({"Br": "Br2", "Cl": "Cl2"}.get(f, f)) for f in (k.split(" [")[0] for k in d)) if d else None  # elemental halogen ref is the X2 solid


def win_text(tp, ion):  # whole phrase, so the prose reads right whether or not the model has a window
    w = WIN.get(tp, {}).get("window")
    if not w: return None, PENDING
    tag = "0 K PBE-D3 grand-potential hull"
    if not w["stable"]: return f"none (at this setting the model lies above the grand-potential hull at every voltage; {tag})", None
    lo = "the metal potential" if w.get("note_red") else f"{w['V_reduction']:.2f} V"  # note_red: no reduction limit above 0 V
    hi = (f"{w['V_oxidation']:.0f} V scan cap, no oxidation limit found below it (not a stability claim)" if w.get("note_ox")
          else f"{w['V_oxidation']:.2f} V")
    red = products(w.get("reduction_products")) or "none above the metal potential"
    ox = products(w.get("oxidation_products")) or "none found"
    return f"down to {lo} and up to {hi} vs {ion}/{ion}^+^ (reduction products: {red}; oxidation products: {ox}) ({tag})", None


def win_num(tp, field):  # short numbers for the abstract: V_reduction, V_oxidation, window_V (width)
    w = WIN.get(tp, {}).get("window")
    if not w: return None, PENDING
    return (w[field], None) if w["stable"] else (None, MISSING)  # no window -> the text using these must change


def gap(tp):
    g = WIN.get(tp, {}).get("pbe_gap")
    return (g["gap_eV"], None) if g else (None, PENDING)


def g4_max_shift():  # largest |tighter-setting margin - base margin| over both models (numeric entries only)
    out = []
    for tp in (NA, LI):
        base = {"dE": C[tp]["dE_vs_listed"], "competing": C[tp]["vs_competing"]}
        for s in C[tp]["tighter_settings"].values():
            out += [abs(v - base[k]) for k, v in s.items() if isinstance(v, (int, float))]
    return (max(out), None) if out else (None, MISSING)


def md_both():  # one MD sentence when both models gave the same result
    (a, ea), (b, eb) = ALIAS["na_md"](), ALIAS["li_md"]()
    if a is None or b is None: return None, ea or eb
    return (f"{a}, in both models" if a == b else f"for TP-NaSbBr4-01, {a}; for TP-LiSbBr4-01, {b}"), ea or eb


ALIAS = {}
for f, tp, ion in (("NaSbBr4", NA, "Na"), ("LiSbBr4", LI, "Li")):  # {{window_NaSbBr4}}, {{gap_NaSbBr4}}, ...
    ALIAS |= {f"window_{f}": lambda tp=tp, ion=ion: win_text(tp, ion), f"gap_{f}": lambda tp=tp: gap(tp),
              f"win_lo_{f}": lambda tp=tp: win_num(tp, "V_reduction"), f"win_hi_{f}": lambda tp=tp: win_num(tp, "V_oxidation"),
              f"win_width_{f}": lambda tp=tp: win_num(tp, "window_V")}
for pre, tp in (("na", NA), ("li", LI)):
    ALIAS |= {f"{pre}_setting": lambda tp=tp: claim(tp, "setting"),
              f"{pre}_listed_id": lambda tp=tp: claim(tp, "best_listed"),
              f"{pre}_dE_listed": lambda tp=tp: claim(tp, "dE_vs_listed"),
              f"{pre}_comp": lambda tp=tp: claim(tp, "vs_competing"),
              f"{pre}_route": lambda tp=tp: (lambda v, e: (sub_formula(v) if v else v, e))(*claim(tp, "competing_route")),
              f"{pre}_dE_listed_k017": lambda tp=tp: claim(tp, "tighter_settings", "80 Ry k0.17", "dE"),
              f"{pre}_comp_k017": lambda tp=tp: claim(tp, "tighter_settings", "80 Ry k0.17", "competing"),
              f"{pre}_dE_listed_100": lambda tp=tp: claim(tp, "tighter_settings", "100 Ry k0.22", "dE"),
              f"{pre}_comp_100": lambda tp=tp: claim(tp, "tighter_settings", "100 Ry k0.22", "competing"),
              f"{pre}_md": lambda tp=tp: (lambda v, e: (v.split(";")[0] if v else v, e))(*claim(tp, "md")),  # no conductivity bound
              f"{pre}_hull_bj": lambda tp=tp: claim(tp, "dispersion", "D3BJ_hull"),
              f"{pre}_hull_none": lambda tp=tp: claim(tp, "dispersion", "noD3_hull"),
              f"{pre}_sbbr3_flip": lambda tp=tp: sbbr3_flip(tp),
              f"{pre}_pending_list": lambda tp=tp: pending_names(tp)}
ALIAS |= {"na_ph_mesh_min": lambda: claim(NA, "phonons", "min_mesh_THz"),
          "na_ph_away_gamma": lambda: claim(NA, "phonons", "min_away_gamma_THz"),
          "na_ph_comm_min": lambda: claim(NA, "phonons", "min_commensurate_THz"),
          "na_ph_tol": lambda: claim(NA, "phonons", "tolerance_THz"),
          "na_ph_supercell": lambda: (lambda v, e: ("x".join(map(str, v)) if v else v, e))(*claim(NA, "phonons", "supercell")),
          "li_phonon_status": li_phonon_status,
          "li_ph_mesh_min": lambda: li_ph("min_mesh_THz"),
          "li_ph_away_gamma": lambda: li_ph("min_away_gamma_THz"),
          "li_ph_comm_min": lambda: li_ph("min_commensurate_THz"),
          "li_ph_supercell": lambda: (lambda v, e: ("x".join(map(str, v)) if v else v, e))(*li_ph("supercell")),
          "li_listed2_gap": li_listed2_gap,
          "na_noD3_hull": lambda: claim(NA, "dispersion", "noD3_hull"),
          "na_hull_ob86": lambda: claim(NA, "dispersion", "OB86", "hull"),  # vdW-DF-optB86b, claims.py nonlocal_check()
          "li_hull_ob86": lambda: claim(LI, "dispersion", "OB86", "hull"),
          "na_dE_listed_ob86": lambda: claim(NA, "dispersion", "OB86", "dE_vs_listed"),
          "li_dE_listed_ob86": lambda: claim(LI, "dispersion", "OB86", "dE_vs_listed"),
          "li_dE_listed_none": lambda: claim(LI, "dispersion", "OB86", "noD3_dE_vs_listed"),  # plain PBE; stored in the OB86 block
          "sbbr3_order_ob86": lambda: resolve("nonlocal.SbBr3_Pnma_minus_R3"),
          "sbbr3_alt_order": lambda: f1("Na-Sb-Br", "agm003218590"),  # SbBr3 cell started from P2_1 2_1 2_1 (expt) minus R-3
          "sbbr3_p1": lambda: f1("Na-Sb-Br", "agm002257875"),
          "n_pending": lambda: n_pending(),
          "n_pending_sbbr3": lambda: n_pending("SbBr3"),
          "n_pending_nabr": lambda: n_pending("NaBr"),
          "g4_max_shift": g4_max_shift,
          "md_both": md_both,
          "pnma_contacts": lambda: (lambda c: (f"{min(c):.2f}–{max(c):.2f}", None) if c else (None, MISSING))(
              BLOCKS["struct"].get("sbbr3_pnma", {}).get("secondary_contacts_3p6_to_4p2_A")),
          "na_tern_min": lambda: na_ternaries(min),
          "na_tern_max": lambda: na_ternaries(max)}


def resolve(key):
    if key in ALIAS: return ALIAS[key]()
    if key.startswith("include:"):
        f = P / "figures/data" / key[8:]
        if not f.exists(): return None, MISSING
        t = f.read_text()
        return (t.split("\n", 1)[1] if t.startswith("# git") else t).strip(), None
    head, _, rest = key.partition(".")
    if head in BLOCKS:
        v = BLOCKS[head]
        for k in rest.split("."): v = v.get(k) if isinstance(v, dict) else None
        return (None, MISSING) if v in (None, "") else (None, PENDING) if str(v).upper().startswith("PENDING") else (v, None)
    if head == "T1":
        row, _, col = rest.partition(".")
        if row not in T1_ROWS or col not in T1_COLS: return None, MISSING
        v = T1[T1_ROWS[row]][T1_COLS[col]].strip()
        return (None, PENDING) if v == "pending" else (v, None)  # "not computed" is a true statement, kept verbatim
    return None, MISSING


def fmt(v, spec):
    if not spec: return str(v)
    if spec.startswith("abs"): v, spec = abs(float(v)), spec[3:]
    try: float(v)
    except (TypeError, ValueError): return str(v)  # non-numeric (e.g. "not computed")
    if not spec: return str(v)
    m = re.search(r"\.(\d+)f$", spec)  # round half away from zero (-4.55 -> -4.6), not float/half-even
    d = Decimal(str(v))
    return format(d.quantize(Decimal(1).scaleb(-int(m.group(1))), rounding=ROUND_HALF_UP) if m else d, spec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draft", action="store_true", help="render anyway; missing keys become [[key]]")
    a = ap.parse_args()
    src = (P / "preprint_v2.md").read_text()
    bad = {}

    def sub(m):
        key, spec = m.group(1).strip(), m.group(2)
        v, err = resolve(key)
        if err: bad[key] = err
        return f"[[{key}]]" if v is None else fmt(v, spec)

    def cond(m):
        key, want, yes, no = m.group(1).strip(), m.group(2), m.group(3), m.group(4)
        v, err = resolve(key)
        if want is not None:
            if err == MISSING: bad[key] = err
            return yes if v is not None and str(v) == want.strip() else no
        if err == PENDING: bad[key] = err  # a pending fallback branch must still block the final render
        empty = v is None or v in ("", "none") or (type(v) is int and v == 0)  # int 0 = a zero count (n_pending); 0.0 eV renders
        return no if empty else yes

    body = r"((?:(?!\{\{\?).)*?)"  # a branch holding no further opener -> innermost conditional first, so nesting works
    pat = re.compile(r"\{\{\?([^}=]+)(?:=([^}]*))?\}\}" + body + r"\{\{:\}\}" + body + r"\{\{/\?\}\}", re.S)
    while (new := pat.sub(cond, src)) != src: src = new
    out = re.sub(r"\{\{([^}|]+)(?:\|([^}]*))?\}\}", sub, src)
    probs = sorted(f"{k}: {e}" for k, e in bad.items())
    if not a.draft:
        if probs: sys.exit("fill_preprint: refusing to write, unresolved keys:\n  " + "\n  ".join(probs))
        (P / "preprint_v2.filled.md").write_text(out); print("wrote paper/preprint_v2.filled.md"); return
    miss = [p for p in probs if p.endswith(MISSING)]
    banner = (f"> DRAFT RENDER (claims.json commit {J.get('commit')}). Not for submission. {len(probs)} unresolved keys "
              f"({len(miss)} missing, shown as [[key]]; the rest PENDING, shown with their current value):\n> "
              + "; ".join(probs) + "\n\n")
    out = re.sub(r"^(---\n.*?\n---\n)", lambda m: m.group(1) + "\n" + banner, out, count=1, flags=re.S)
    (P / "preprint_v2.draft.md").write_text(out)
    print(f"wrote paper/preprint_v2.draft.md ({len(probs)} unresolved)\n  " + "\n  ".join(probs))


if __name__ == "__main__":
    main()
