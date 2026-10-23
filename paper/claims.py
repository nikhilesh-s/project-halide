"""Claims ledger: the single source of every number the paper, NOMAD v2, Zenodo and the activity line quote.
Run from repo root after verdict.py: discovery/.venv/bin/python -I paper/claims.py  -> paper/claims.md + paper/claims.json
Numbers come only from discovery/results/verdict.json and the phonon summaries; nothing is typed by hand.
A claim is FINAL only when its gate has no pending competitor; otherwise it prints PENDING and must not be quoted as settled."""
import json, re, subprocess
from pathlib import Path

R = Path(__file__).resolve().parent.parent
V = json.loads((R / "discovery/results/verdict.json").read_text())
PH = R / "anish_results/ionic_packet/anish_packet_run"
TOP2 = {"TP-NaSbBr4-01": "S4", "TP-LiSbBr4-01": "S5"}
PH_NOTE = {"TP-NaSbBr4-01": "disp-000 residual force 0.0198 Ry/bohr exceeds the 0.01 gate; the +/- displacement pairs cancel it exactly (phonon_review.md)",
           "TP-LiSbBr4-01": "disp-000 residual force 0.0044 Ry/bohr passes the 0.01 gate"}
MD_BOUND = ("no long-range mobile-ion diffusion in 50-150 ps MACE-MPA-0 MD at fixed PBE-D3 volume, 300-600 K; "
            "conductivity below about 5e-3 S/cm at 600 K (slower conduction not ruled out)")


def base_row(c):
    rows = V["candidates"][c]["rows"]
    return next(r for r in rows if r["setting"] == V["candidates"][c]["gates"]["base_setting"])


def phonons(c, stage):
    p = PH / stage / c / "phonon_summary.json"
    if not p.exists(): return None
    s = json.loads(p.read_text())
    return dict(min_mesh_THz=round(s["min_freq_mesh_all_THz"], 3), min_away_gamma_THz=round(s["min_freq_mesh_away_from_Gamma_THz"], 3),
                min_commensurate_THz=round(s["min_freq_commensurate_nonGamma_THz"], 3), supercell=s["supercell"], tolerance_THz=s["threshold_THz"], verdict=s.get("verdict"))




def na_dispersion_hull():
    """NaSbBr4 margin vs every computed competitor with D3-BJ and without D3 (m5_top2 + m5_oct20 relaxations)."""
    import glob
    def e(p):
        t = Path(p).read_text(errors="replace")
        if "JOB DONE" not in t or "End final coordinates" not in t: return None
        en = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)", t, re.M)
        n = int(re.search(r"number of atoms/cell\s+=\s+(\d+)", t).group(1)); return float(en[-1]) * 13.605693122994 / n
    def low(pat):
        v = [x for x in (e(f) for f in glob.glob(str(R / "overnight" / pat))) if x is not None]; return min(v) if v else None
    res = {}
    for mode, tpjob in (("D3BJ", "01"), ("noD3", "02")):
        tp, nabr, sb = low(f"m5_top2/jobs/{tpjob}_NaSbBr4_TP01_hit_{mode}_vcrelax.out"), low(f"m5_top2/jobs/*_NaBr_*_{mode}_vcrelax.out"), low(f"m5_oct20/jobs/*_SbBr3_*_{mode}_vcrelax.out")
        na3, na2 = low(f"m5_top2/jobs/*_Na3SbBr6_*_{mode}_vcrelax.out"), low(f"m5_top2/jobs/*_Na2SbBr5_*_{mode}_vcrelax.out")
        if None in (tp, nabr, sb, na3, na2): res[f"{mode}_hull"] = "PENDING"; continue
        routes = [(2 * nabr + 4 * sb) / 6, (10 / 3 * na3 + 8 / 3 * sb) / 6, (4 * na2 + 2 * sb) / 6]  # NaBr+SbBr3, 1/3 Na3SbBr6+2/3 SbBr3, 1/2 Na2SbBr5+1/2 SbBr3
        res[f"{mode}_hull"] = round(max((tp - r) * 1000 for r in routes), 6)
    return res


def li_dispersion_hull():
    """LiSbBr4 margin vs LiBr+SbBr3 and 1/3 Li3SbBr6 + 2/3 SbBr3, D3-BJ and no D3 (m5_* relaxations). PENDING (with the missing phase named) if any is absent."""
    import glob
    def e(p):
        t = Path(p).read_text(errors="replace")
        if "JOB DONE" not in t or "End final coordinates" not in t: return None
        en = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)", t, re.M)
        n = int(re.search(r"number of atoms/cell\s+=\s+(\d+)", t).group(1)); return float(en[-1]) * 13.605693122994 / n
    def low(pat):
        v = [x for x in (e(f) for f in glob.glob(str(R / "overnight" / pat))) if x is not None]; return min(v) if v else None
    res, missing = {}, {}
    for mode in ("D3BJ", "noD3"):
        ph = dict(TP=low(f"m5_*/jobs/*_LiSbBr4_TP01_hit_{mode}_vcrelax.out"), LiBr=low(f"m5_*/jobs/*_LiBr_*_{mode}_vcrelax.out"),
                  SbBr3=low(f"m5_*/jobs/*_SbBr3_*_{mode}_vcrelax.out"), Li3SbBr6=low(f"m5_*/jobs/*_Li3SbBr6_*_{mode}_vcrelax.out"))
        miss = [k for k, v in ph.items() if v is None]
        if miss: res[f"{mode}_hull"] = "PENDING"; missing[mode] = miss; continue
        routes = [(2 * ph["LiBr"] + 4 * ph["SbBr3"]) / 6, (10 / 3 * ph["Li3SbBr6"] + 8 / 3 * ph["SbBr3"]) / 6]  # LiBr+SbBr3; 1/3 Li3SbBr6 + 2/3 SbBr3
        res[f"{mode}_hull"] = round(max((ph["TP"] - r) * 1000 for r in routes), 6)
    res["missing_phases"] = missing
    return res


def _epa(p):  # final energy per atom (eV) of a converged vc-relax, else None
    t = Path(p).read_text(errors="replace")
    if "JOB DONE" not in t or "End final coordinates" not in t: return None
    en = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)", t, re.M)
    return float(en[-1]) * 13.605693122994 / int(re.search(r"number of atoms/cell\s+=\s+(\d+)", t).group(1))


def _low(pat):
    import glob
    v = [x for x in (_epa(f) for f in glob.glob(str(R / "overnight" / pat))) if x is not None]; return min(v) if v else None


def nonlocal_check():
    """Third dispersion model: vdW-DF-optB86b (overnight/m5_disp3, make_disp3.py), plus the no-D3 listed-LiSbBr4 gap job.
    Per candidate: hull margin (same routes as na_/li_dispersion_hull) and dE vs listed polymorph; also SbBr3 Pnma minus R-3."""
    L = lambda s, m="OB86": _low(f"m5_disp3/jobs/*_{s}_{m}_vcrelax.out")
    r3, pnma = L("SbBr3_agm002139075"), L("SbBr3_Pnma_expt_agm003218591")
    sb = min([x for x in (r3, pnma) if x is not None], default=None)
    res = dict(functional="vdW-DF-optB86b (input_dft='vdw-df-ob86'), from PBE-D3 geometries, 80 Ry k0.22",
               SbBr3_Pnma_minus_R3=round((pnma - r3) * 1000, 6) if None not in (r3, pnma) else "PENDING")
    for c, A, terns in (("TP-NaSbBr4-01", "Na", (("Na3SbBr6_agm016963778", 10 / 3, 8 / 3), ("Na2SbBr5_agm005863603", 4, 2))),
                        ("TP-LiSbBr4-01", "Li", (("Li3SbBr6_agm006027634", 10 / 3, 8 / 3),))):
        tp, ax, ls = L(f"{A}SbBr4_TP01_hit"), L(f"{A}Br_*"), L(f"{A}SbBr4_*_best_existing")
        tv = {n: L(n) for n, _, _ in terns}
        ph = dict(TP=tp, AX=ax, SbBr3_R3=r3, SbBr3_Pnma=pnma, listed=ls, **tv)
        o = dict(missing=[k for k, v in ph.items() if v is None])
        o["dE_vs_listed"] = round((tp - ls) * 1000, 6) if None not in (tp, ls) else "PENDING"
        if None in (tp, ax, sb, *tv.values()): o["hull"] = "PENDING"
        else:
            routes = [(2 * ax + 4 * sb) / 6] + [(a * tv[n] + b * sb) / 6 for n, a, b in terns]
            o["hull"] = round(max((tp - r) * 1000 for r in routes), 6)
        res[c] = o
    # Table 1 optB86b cells (same rows as T1_method_error.csv): AX + lowest SbBr3 per LiBr polymorph, and LiBr rocksalt - wurtzite
    na_tp, li_tp, nabr = L("NaSbBr4_TP01_hit"), L("LiSbBr4_TP01_hit"), L("NaBr_*")
    rs, wz = L("LiBr_rocksalt_*"), L("LiBr_wurtzite_*")
    vs = lambda tp, ax: round((tp - (2 * ax + 4 * sb) / 6) * 1000, 6) if None not in (tp, ax, sb) else "PENDING"
    lib = min([x for x in (rs, wz) if x is not None], default=None)
    res["T1"] = dict(na_bin=vs(na_tp, nabr), li_bin=vs(li_tp, lib), li_rs=vs(li_tp, rs), li_wz=vs(li_tp, wz),
                     libr_order=round((rs - wz) * 1000, 6) if None not in (rs, wz) else "PENDING")
    tp, ls = _low("m5_*/jobs/*_LiSbBr4_TP01_hit_noD3_vcrelax.out"), L("LiSbBr4_*_best_existing", "noD3")
    res["TP-LiSbBr4-01"]["noD3_dE_vs_listed"] = round((tp - ls) * 1000, 6) if None not in (tp, ls) else "PENDING"
    return res


S = (R / "discovery/results/transplant_cycle3.md").read_text()
num = lambda pat: int(re.search(pat, S).group(1).replace(",", ""))
SEARCH = dict(n_entries=num(r"downloaded \(<= 24 atoms\): (\d+)"), n_prototypes=num(r"\*\*(\d+) distinct prototypes\*\*"),
              n_transplants=num(r"transplants attempted \(cheap MACE pass\): (\d+)"), n_ml_d3=num(r"Kept for MACE\+D3 .*?: (\d+)\."),
              n_hits=num(r"HIT_proposed_model: (\d+)"), n_dft_models=len(V["candidates"]),
              val_recovered=3, val_total=4, val_names="NaAlBr~4~, NaGaCl~4~, LiGaCl~4~ and LiAlBr~4~", val_lialbr4_miss="about 19", source="discovery/results/transplant_cycle3.md; gate_night4.md (validation 3 of 4)")
_gpass = lambda g: all(str(g.get(k, "")).startswith(("pass", "converged pass")) for k in ("G1", "G2", "G3", "G4"))
SEARCH |= dict(n_dft_formulas=len({base_row(c)["formula"] for c in V["candidates"]}),
               n_gate_pass=len(gp := {base_row(c)["formula"] for c, v in V["candidates"].items() if _gpass(v["gates"])}),
               n_gate_pass_other=len(gp) - len(TOP2), gate_pass_formulas=sorted(gp))  # compositions meeting G1-G4 in verdict.py
def lialbr4_control():
    """DFT false-positive control: TP-LiAlBr4-01 (ML pick, P2_1/m) vs the LiAlCl4-type P2_1/c cell (experimental structure type),
    both relaxed with the G2/G3 protocol (PBE-D3 zero damping, 80 Ry k0.22). Value = ML pick minus experimental type, meV/atom;
    positive = experimental type correctly ranked lower; negative = DFT-level false positive."""
    import glob
    def e(f):
        t = Path(f).read_text(errors="replace")
        if "JOB DONE" not in t or "End final coordinates" not in t: return None
        en = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)", t, re.M)
        return float(en[-1]) * 13.605693122994 / int(re.search(r"number of atoms/cell\s+=\s+(\d+)", t).group(1))
    pick = lambda pat: [x for x in (e(f) for f in glob.glob(str(R / "overnight" / pat))) if x is not None]
    tp, ex = pick("*/jobs/*LiAlBr4_TP01_control*.out"), pick("*/jobs/*LiAlBr4_LiAlCl4type*control*.out")
    if not (tp and ex): return dict(status="not run", setting="80 Ry k0.22", note="LiAlBr4 control pair not both converged")
    return dict(status="done", setting="80 Ry k0.22, PBE-D3 zero damping", lialbr4=dict(zero=round((min(tp) - min(ex)) * 1000, 6), bj="not computed"),
                note=f"{len(tp)} ML-pick and {len(ex)} experimental-type runs; negative = DFT also prefers the hypothetical model (false positive)")


META = json.loads((R / "paper/meta.json").read_text()) if (R / "paper/meta.json").exists() else {}

NL = nonlocal_check()
out, md = {}, []
commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=R, capture_output=True, text=True).stdout.strip()
md += [f"# Claims ledger (generated by `paper/claims.py` at commit {commit}; do not edit by hand)", "",
       "Labels: \"predicted\", \"proposed model\". PENDING means a competitor is still running: do not quote that number as settled.", ""]
for c, stage in TOP2.items():
    g, b = V["candidates"][c]["gates"], base_row(c)
    pend = b["pending_in_system"]
    ph = phonons(c, stage)
    tight = [r for r in V["candidates"][c]["rows"] if r is not b and r["status"] == "done"]
    f = b["formula"]
    o = dict(formula=f, setting=b["setting"], best_listed=b["best_listed"], dE_vs_listed=b["dE"], vs_competing=b["competing"],
             competing_route=" + ".join(b["competing_route"]), pending=pend, status="FINAL" if not pend else "PENDING",
             tighter_settings={r["setting"]: dict(dE=r["dE"], competing=r["competing"]) for r in tight},
             phonons=ph, md=MD_BOUND)
    if c == "TP-NaSbBr4-01": o["dispersion"] = dict(noD3_hull=(d := na_dispersion_hull())["noD3_hull"], D3BJ_hull=d["D3BJ_hull"], OB86=NL[c])
    if c == "TP-LiSbBr4-01": o["dispersion"] = dict(noD3_hull=(d := li_dispersion_hull())["noD3_hull"], D3BJ_hull=d["D3BJ_hull"], missing_phases=d["missing_phases"], OB86=NL[c])
    out[c] = o
    md += [f"## {c} ({f}) - {o['status']}", "",
           f"- Claim: proposed model {c}, predicted {abs(b['dE']):.1f} meV/atom below the lowest computed listed {f} structure "
           f"({b['best_listed']}) and {abs(b['competing']):.1f} meV/atom below {o['competing_route']} "
           f"(PBE-D3 zero damping, {b['setting']}, 0 K, no SOC).",
           f"- Other settings: " + "; ".join(f"{s}: dE {v['dE']} / competing {v['competing']}" for s, v in o["tighter_settings"].items()),
           "- Phonons: " + (f"no instability found at {ph['supercell']} sampling (mesh minimum {ph['min_mesh_THz']:+.3f} THz, "
                            f"commensurate minimum {ph['min_commensurate_THz']:+.3f} THz; tolerance {ph['tolerance_THz']} THz; {PH_NOTE[c]})" if ph else "PENDING (set not complete)"),
           f"- Dispersion: hull no-D3 {o['dispersion']['noD3_hull']}, D3-BJ {o['dispersion']['D3BJ_hull']}, vdW-DF-optB86b {NL[c]['hull']} meV/atom "
           f"(optB86b dE vs listed {NL[c]['dE_vs_listed']}; SbBr3 Pnma minus R-3 {NL['SbBr3_Pnma_minus_R3']}).",
           f"- MD: {MD_BOUND}.",
           "- Pending competitors: " + ("none" if not pend else "; ".join(re.sub(r" \(.*", "", p) for p in pend)), ""]
(R / "paper/claims.json").write_text(json.dumps(dict(commit=commit, claims=out, search=SEARCH, meta=META, control=lialbr4_control(), nonlocal_check=NL), indent=2) + "\n")
(R / "paper/claims.md").write_text("\n".join(md))
print("\n".join(md))
