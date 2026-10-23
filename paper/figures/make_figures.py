"""Paper figures F1-F3 and table T1 for Project Halide (TP-NaSbBr4-01, TP-LiSbBr4-01: PROPOSED MODELS, predicted).

Run from anywhere:  discovery/.venv/bin/python -I paper/figures/make_figures.py
Everything is rebuilt from repo data; nothing outside paper/figures/ is written. Energies/phases come from
discovery/scripts/verdict.py (parse, analyse, decompose are imported, not re-implemented). Every output carries the git commit
and the mtime of discovery/results/verdict.json.
"""
import glob, importlib.util, subprocess, sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core import Composition
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SA

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "discovery/scripts"))
import verdict as V  # noqa: E402  (parse, analyse, decompose, mi.final_structure, RY)

DATA = HERE / "data"; DATA.mkdir(exist_ok=True)
BASE = "80 Ry k0.22"                    # setting every TP/listed/binary comparison in the paper is made at (verdict.py base setting)
RUN = ROOT / "anish_results/ionic_packet/anish_packet_run"
PACKET = ROOT / "anish_packet"
TPS = {"Na": dict(tp="TP-NaSbBr4-01", f="NaSbBr4", system="Br-Na-Sb", title="Na-Sb-Br", S="S4"),
       "Li": dict(tp="TP-LiSbBr4-01", f="LiSbBr4", system="Br-Li-Sb", title="Li-Sb-Br", S="S5")}
METHODS = {"zero": "PBE-D3 (zero damping)", "bj": "PBE-D3(BJ)", "none": "PBE (no D3)"}
# Okabe-Ito (colorblind safe)
C_TP, C_LISTED, C_OTHER, C_GREY, C_BAND = "#D55E00", "#0072B2", "#56B4E9", "#7F7F7F", "#DDDDDD"

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9, "legend.fontsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "savefig.dpi": 300, "figure.dpi": 100, "axes.linewidth": 0.6, "lines.linewidth": 0.9})


def git(*a):
    return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True).stdout.strip()


COMMIT = git("describe", "--always", "--dirty") or "no-git"
VJ = ROOT / "discovery/results/verdict.json"
VTIME = datetime.fromtimestamp(VJ.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
STAMP = f"git {COMMIT} | verdict.json {VTIME}"


def save(fig, name, dy=-0.03):
    fig.savefig(HERE / f"{name}.png", metadata={"Description": STAMP}, bbox_inches="tight")
    fig.savefig(HERE / f"{name}.pdf", metadata={"Subject": STAMP}, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name, ".png/.pdf")


# ---------------------------------------------------------------- data (verdict.py parse, all dispersion variants)
def load_rows():
    parked = [f for f in glob.glob(str(ROOT / "overnight/parked/*/*.in")) if Path(f).with_suffix(".out").exists()]
    rows = []
    for f in sorted(glob.glob(str(ROOT / "overnight/*/jobs/*.in")) + parked):
        if Path(f).parent.parent.name in V.SKIP: continue
        t = Path(f).read_text(); r = V.parse(f); r["path"] = Path(f).relative_to(ROOT).as_posix()
        r["method"] = "ob86" if "ob86" in t.lower() else (("bj" if "dftd3_version" in t else "zero") if "vdw_corr='dft-d3'" in t else "none")   # optB86b jobs have no vdw_corr: not "no D3"
        rows.append(r)
    return rows


ALL = load_rows()
PATHS = {r["job"]: r["path"] for r in ALL}   # job id (as verdict.py names it) -> repo-relative .in path
ROWS = {m: [r for r in ALL if r["method"] == m] for m in METHODS}
BEST, CANDS = {}, {}
for m in METHODS:
    BEST[m], CANDS[m] = V.analyse(ROWS[m])


def tp_row(m, key):
    return next((v for v in CANDS[m].get(TPS[key]["tp"], []) if v["setting"] == BASE), None)


def lowest(m, source=None, formula=None, tp=False):
    """Lowest-energy finished row at BASE for a source id or formula (per-atom energy)."""
    c = [r for r in ROWS[m] if r["status"] == "done" and r["setting"] == BASE and r["tp"] == tp
         and (source is None or r["source"] == source) and (formula is None or r["formula"] == formula)]
    return min(c, key=lambda r: r["e_ry"] / r["natoms"]) if c else None


# consistency with the committed verdict.json (D3 zero): warn, do not stop
vj = __import__("json").loads(VJ.read_text())
for kk, k in TPS.items():
    mine, theirs = tp_row("zero", kk), next(
        (v for v in vj["candidates"][k["tp"]]["rows"] if v["setting"] == BASE), None)
    for key in ("dE", "binaries", "competing"):
        if mine and theirs and isinstance(mine[key], float) and isinstance(theirs[key], float) and abs(mine[key] - theirs[key]) > 0.011:
            print(f"WARNING: {k['tp']} {key}: fresh parse {mine[key]} != verdict.json {theirs[key]} (verdict.json older than job outputs?)")


# cells whose relaxed symmetry differs from the cell they were started from (relaxed cell is P2_12_12_1 at symprec 0.001, Pnma at 0.05)
START_NOTE = {"agm003218590": "  (P2$_1$2$_1$2$_1$ start; Pnma at symprec 0.05)"}


# ---------------------------------------------------------------- F1 competing-phase energies
def f1():
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 8.4), gridspec_kw=dict(hspace=0.28))
    pnma, r3 = lowest("zero", "agm003218591"), lowest("zero", "agm002139075")
    err = (pnma["e_ry"] / pnma["natoms"] - r3["e_ry"] / r3["natoms"]) * V.RY if pnma and r3 else None
    export = []
    for ax, key in zip(axes, TPS):
        t = TPS[key]; els = set(t["system"].split("-"))
        rows = [r for r in ROWS["zero"] if r["setting"] == BASE and {e.symbol for e in Composition(r["formula"])} <= els]
        refs = {}   # lowest computed binary per formula (Ry/atom)
        for r in rows:
            if r["status"] == "done" and not r["tp"] and len(Composition(r["formula"])) == 2:
                e = r["e_ry"] / r["natoms"]; refs[r["formula"]] = min(e, refs.get(r["formula"], 9e9))
        listed = lowest("zero", formula=t["f"]); items = {}
        for r in rows:
            k = (r["formula"], r["source"])
            if r["status"] == "done":
                d = V.decompose(r["formula"], refs) if len(Composition(r["formula"])) > 1 else None
                if d is None: continue
                de = (r["e_ry"] / r["natoms"] - d[0]) * V.RY
                if k not in items or items[k].get("de") is None or de < items[k]["de"]:
                    items[k] = dict(de=de, tp=r["tp"], sg=r["sg"], job=r["job"], listed=bool(listed and r["job"] == listed["job"]))
            elif k not in items:
                items[k] = dict(de=None, tp=r["tp"], sg=None, job=r["job"], listed=False)
        order = sorted(items.items(), key=lambda kv: (kv[1]["de"] is None, kv[1]["de"] if kv[1]["de"] is not None else 0))
        n = len(order); xs = [v["de"] for _, v in order if v["de"] is not None]
        lo, hi = min(xs + [-err if err else 0]) - 6, max(xs + [err if err else 0]) + 22
        ax.axvspan(-err, err, color=C_BAND, zorder=0, lw=0)
        ax.axvline(0, color="k", lw=0.6, zorder=1)
        for i, ((f, src), v) in enumerate(order):
            y = n - 1 - i
            if v["de"] is None:
                ax.plot(0.955, y, "o", mfc="none", mec=C_GREY, ms=5.5, mew=1, transform=ax.get_yaxis_transform(), clip_on=False)
                ax.text(0.945, y, "pending", ha="right", va="center", fontsize=6.5, color=C_GREY, transform=ax.get_yaxis_transform())
                lab_c = C_GREY
            else:
                c = C_TP if v["tp"] else C_LISTED if v["listed"] else C_OTHER if abs(v["de"]) > 1e-9 else C_GREY
                ax.plot(v["de"], y, "o", color=c, ms=6 if v["tp"] else 5, zorder=3, mec="k" if v["tp"] else "none", mew=0.6)
                ax.text(v["de"] + 1.2, y, ("0" if abs(v["de"]) < 1e-9 else f"{v['de']:+.1f}"), va="center", fontsize=6.5, color="k", zorder=4)
                lab_c = "k"
            ax.text(-0.012, y, f"{f}  {src if v['tp'] else src.replace('agm', 'agm')}" + (START_NOTE.get(src) or (f"  ({v['sg']})" if v["sg"] else "")),
                    ha="right", va="center", fontsize=6.8, color=lab_c, fontweight="bold" if v["tp"] else "normal", transform=ax.get_yaxis_transform())
            export.append(dict(system=t["title"], formula=f, source=src, space_group=("Pnma (P212121 start; symprec 0.05)" if src in START_NOTE else v["sg"] or ""), job=v["job"], dE_vs_binaries_meV_atom="" if v["de"] is None else f"{v['de']:.6f}",
                               status="pending" if v["de"] is None else "done", role="TP" if v["tp"] else "lowest listed" if v["listed"] else ""))
        ax.set_ylim(-0.7, n - 0.3); ax.set_xlim(lo, hi); ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        bins = [b.replace("3", "$_3$") for b in sorted(refs)]
        ax.set_xlabel(f"Energy relative to lowest computed binaries ({' + '.join(bins)}), meV/atom")
        ax.set_title(f"{t['title']}: {t['tp']} (proposed model)", loc="left", fontweight="bold")
        if err: ax.text(err + 0.6, n - 0.45, f"$\\pm${err:.1f}: SbBr$_3$ Pnma vs R-3\nPBE-D3 ordering error", fontsize=6, color="#555", va="top")
    h = [Line2D([], [], marker="o", ls="", color=C_TP, mec="k", ms=6, label="proposed model (TP)"),
         Line2D([], [], marker="o", ls="", color=C_LISTED, ms=5, label="lowest computed listed polymorph, same composition"),
         Line2D([], [], marker="o", ls="", color=C_OTHER, ms=5, label="other computed phase"),
         Line2D([], [], marker="o", ls="", color=C_GREY, ms=5, label="binary reference (0)"),
         Line2D([], [], marker="o", ls="", mfc="none", mec=C_GREY, ms=5.5, label="pending (no energy yet)")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.62, -0.035))
    fig.text(0.5, -0.075, f"PBE-D3 (zero damping), QE, {BASE.replace('k', 'k-spacing ')}, 0 K, no SOC. Negative = below the binaries. Predicted, not experimentally verified.",
             ha="center", va="top", fontsize=6.5, color="#444")
    fig.subplots_adjust(left=0.40, right=0.97, top=0.96, bottom=0.05)
    save(fig, "F1_competing_phases")
    import csv
    with open(DATA / "F1_competing_phases.csv", "w", newline="") as fh:
        fh.write(f"# {STAMP}; setting {BASE}; PBE-D3 zero damping; energies relative to the lowest computed binaries\n")
        w = csv.DictWriter(fh, fieldnames=list(export[0])); w.writeheader(); w.writerows(export)


# ---------------------------------------------------------------- F2 phonons
PATH = [[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [0, .5, 0], [0, 0, 0], [0, 0, .5], [.5, 0, .5], [.5, .5, .5], [0, .5, .5], [0, 0, 0]]  # same generic path as post_phonons.py
PLAB = ["$\\Gamma$", "½00", "½½0", "0½0", "$\\Gamma$", "00½", "½0½", "½½½", "0½½", "$\\Gamma$"]


def load_phonopy(key):
    """Returns (phonopy object, note) or (None, reason). Li fills in automatically once all displacement runs finished."""
    import phonopy
    t = TPS[key]; d = RUN / t["S"] / t["tp"]
    outs = sorted(d.glob("disp-*.out")); done = [o for o in outs if "JOB DONE" in o.read_text(errors="replace")]
    y = d / "phonopy_params.yaml"
    if y.exists(): return phonopy.load(y, log_level=0), d, len(done), len(outs)
    dis = PACKET / t["S"] / t["tp"] / "phonopy_disp.yaml"
    if dis.exists() and outs and len(done) == len(outs) and len(outs) == 19:       # 19/19: same recipe as post_phonons.py, in memory only
        spec = importlib.util.spec_from_file_location("post_phonons", PACKET / "S1/post_phonons.py"); pp = importlib.util.module_from_spec(spec); spec.loader.exec_module(pp)
        ph = phonopy.load(dis, produce_fc=False, log_level=0); nat = len(ph.supercell); f0, _ = pp.qe_forces(d / "disp-000.out", nat)
        ph.forces = np.array([pp.qe_forces(d / f"disp-{i:03d}.out", nat)[0] - f0 for i in range(1, len(ph.displacements) + 1)])
        ph.produce_force_constants(); ph.symmetrize_force_constants()
        return ph, d, len(done), len(outs)
    return None, d, len(done), len(outs)


def total_force(d):
    import re
    m = re.findall(r"Total force =\s+([\d.]+)", (d / "disp-000.out").read_text(errors="replace"))
    return float(m[-1]) if m else None


def f2():
    fig = plt.figure(figsize=(7.0, 6.6)); gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], hspace=0.38, wspace=0.05)
    for row, key in enumerate(TPS):
        t = TPS[key]; ax, axd = fig.add_subplot(gs[row, 0]), fig.add_subplot(gs[row, 1])
        ph, d, ndone, ntot = load_phonopy(key)
        if ph is None:
            for a in (ax, axd): a.set_xticks([]); a.set_yticks([]); [s.set_visible(True) for s in a.spines.values()]
            ax.text(0.5, 0.5, f"PLACEHOLDER: {t['tp']} phonons not yet available\n{ndone}/{ntot or 19} displacement runs finished (S5, Anish).\n"
                    "Rerun make_figures.py once 19/19 are done;\nthe panel then fills automatically.", ha="center", va="center", transform=ax.transAxes, color=C_TP, fontsize=8)
            ax.set_title(f"{t['tp']}: pending", loc="left", fontweight="bold"); axd.axis("off"); continue
        ph.run_band_structure([np.linspace(a, b, 41) for a, b in zip(PATH[:-1], PATH[1:])])
        bs = ph.get_band_structure_dict(); xt = [0]
        for dist, fr in zip(bs["distances"], bs["frequencies"]):
            ax.plot(dist, fr, color="#0072B2", lw=0.8)
            xt.append(dist[-1])
        ph.run_mesh([16, 16, 16]); mf = ph.get_mesh_dict()["frequencies"]; ph.run_total_dos(); dos = ph.get_total_dos_dict()
        fmax = max(mf.max(), max(np.max(f) for f in bs["frequencies"]))
        for a in (ax, axd): a.axhline(-0.3, color=C_TP, ls="--", lw=0.9); a.set_ylim(-1.7, fmax * 1.04)
        ax.axhline(0, color="k", lw=0.4)
        for xv in xt: ax.axvline(xv, color=C_BAND, lw=0.5, zorder=0)
        ax.set_xticks(xt); ax.set_xticklabels(PLAB); ax.set_xlim(0, xt[-1]); ax.set_ylabel("Frequency (THz)")
        axd.fill_betweenx(dos["frequency_points"], 0, dos["total_dos"], color=C_OTHER, lw=0); axd.plot(dos["total_dos"], dos["frequency_points"], color="#0072B2", lw=0.6)
        axd.axhline(0, color="k", lw=0.4); axd.set_yticklabels([]); axd.set_xticks([]); axd.set_xlabel("DOS"); axd.set_xlim(left=0)
        ax.text(0.99, -0.34, "-0.3 THz tolerance", color=C_TP, fontsize=6.5, ha="right", va="top", transform=ax.get_yaxis_transform())
        tf = total_force(d); sc = np.diag(ph.supercell_matrix).astype(int)
        ax.set_title(f"{t['tp']}: {len(ph.supercell)}-atom {sc[0]}x{sc[1]}x{sc[2]} supercell, mesh minimum {mf.min():+.3f} THz", loc="left", fontweight="bold", fontsize=8)
        ax.text(0.01, 0.02, f"disp-000 total force {tf:.4f} Ry/bohr (gate 0.01: {'met' if tf < 0.01 else 'NOT met'}); disp-000 forces subtracted; no LO-TO" if tf else "",
                transform=ax.transAxes, fontsize=6, color="#444", va="bottom")
    fig.text(0.5, 0.01, "Finite-displacement phonopy, PBE-D3 forces (QE 80 Ry). Generic path through zone-boundary points (not a symmetry-standard path). Predicted, 0 K harmonic.",
             ha="center", va="top", fontsize=6.3, color="#444")
    save(fig, "F2_phonons")


# ---------------------------------------------------------------- F3 powder XRD
def structure_of(row):
    inp = ROOT / PATHS[row["job"]]
    return V.mi.final_structure(inp.with_suffix(".out") if row["calc"] == "vc-relax" else inp)


def f3():
    import csv
    calc = XRDCalculator(wavelength="CuKa"); fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True, gridspec_kw=dict(wspace=0.06))
    for ax, key in zip(axes, TPS):
        t = TPS[key]; tpr = lowest("zero", formula=t["f"], tp=True); lst = lowest("zero", formula=t["f"])
        pats = []
        for label, r, sign, col in ((t["tp"], tpr, 1, C_TP), (f"{lst['source']} (lowest listed)", lst, -1, C_LISTED)):
            s = structure_of(r); p = calc.get_pattern(s, scaled=True, two_theta_range=(5, 60))
            ax.vlines(p.x, 0, sign * np.asarray(p.y), color=col, lw=1.0, label=label)
            pats.append((label, r, s, p))
        ax.axhline(0, color="k", lw=0.5); ax.set_xlim(5, 60); ax.set_ylim(-112, 112); ax.set_xlabel("2$\\theta$ (deg), Cu K$\\alpha$")
        ax.set_title(f"{t['f']}", loc="left", fontweight="bold"); ax.legend(frameon=False, loc="upper right", fontsize=6.5)
        ax.set_yticks([-100, -50, 0, 50, 100]); ax.set_yticklabels(["100", "50", "0", "50", "100"])
        with open(DATA / f"F3_xrd_{t['f']}.csv", "w", newline="") as fh:
            fh.write(f"# {STAMP}; simulated Cu Ka (pymatgen XRDCalculator), 5-60 deg 2theta, intensity scaled to 100; relaxed {BASE} PBE-D3 cells\n")
            w = csv.writer(fh); w.writerow(["structure", "job", "space_group", "two_theta_deg", "d_A", "rel_intensity", "hkl"])
            for label, r, s, p in pats:
                for x, y, dd, h in zip(p.x, p.y, p.d_hkls, p.hkls):
                    w.writerow([label, r["job"], r["sg"], f"{x:.4f}", f"{dd:.4f}", f"{y:.3f}", str(tuple(h[0]["hkl"]))])
    axes[0].set_ylabel("Simulated intensity (%)")
    fig.text(0.5, -0.02, "Simulated, not measured. Proposed model up, lowest computed listed polymorph mirrored down.", ha="center", va="top", fontsize=6.5, color="#444")
    save(fig, "F3_xrd", dy=-0.07)


# ---------------------------------------------------------------- F5 simulated XRD of the two predicted structures (stacked, hkl-labelled)
def f5():
    import csv
    from pymatgen.core import Structure
    calc = XRDCalculator(wavelength="CuKa"); fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True, gridspec_kw=dict(hspace=0.3))
    rows = []
    for ax, key in zip(axes, TPS):
        t = TPS[key]; r = lowest("zero", formula=t["f"], tp=True); s = structure_of(r)
        cif = Structure.from_file(HERE / f"cif/{t['tp']}_{t['f']}_PBE-D3_relaxed.cif")   # shipped CIF must be the final vc-relax geometry
        assert StructureMatcher(ltol=0.01, stol=0.02, angle_tol=0.5, scale=False).fit(cif, s), t["tp"]
        p = calc.get_pattern(s, scaled=True, two_theta_range=(10, 60))
        ax.vlines(p.x, 0, p.y, color=C_TP, lw=1.2); ax.set_ylim(0, 125); ax.set_xlim(10, 60)
        # pymatgen leaves near-coincident reflections separate: merge within 0.2 deg, label the strongest member of the 6 strongest groups
        g = []
        for i in np.argsort(p.x):
            if g and p.x[i] - p.x[g[-1][0]] < 0.2: g[-1].append(i)
            else: g.append([i])
        top = sorted(g, key=lambda q: -sum(p.y[i] for i in q))[:6]
        i0 = max(top[0], key=lambda i: p.y[i]); h0 = p.hkls[i0][0]["hkl"]   # check: strongest peak is Bragg's law on a lattice d_hkl
        assert abs(p.d_hkls[i0] - s.lattice.d_hkl(h0)) < 1e-6 and abs(p.x[i0] - 2 * np.degrees(np.arcsin(calc.wavelength / 2 / s.lattice.d_hkl(h0)))) < 0.02
        done = []
        for q in top:
            i = max(q, key=lambda j: p.y[j]); hk = "".join(f"{abs(x)}" if x >= 0 else f"$\\bar{{{-x}}}$" for x in p.hkls[i][0]["hkl"])
            ax.annotate(f"({hk})", (p.x[i], p.y[i]), xytext=(0, 2 + 8 * sum(abs(p.x[i] - x) < 2 for x in done)), textcoords="offset points", ha="center", fontsize=6.5); done.append(p.x[i])
        ax.set_title(f"{t['tp']}  {t['f']} ({r['sg']})", loc="left", fontweight="bold"); ax.set_ylabel("Simulated intensity (%)")
        rows += [[f"{x:.4f}", f"{d:.4f}", str(tuple(h[0]["hkl"])), f"{y:.3f}", t["tp"]] for x, y, d, h in zip(p.x, p.y, p.d_hkls, p.hkls)]
    axes[1].set_xlabel("2$\\theta$ (deg), Cu K$\\alpha$ (1.5418 Å)")
    fig.text(0.5, 0.0, "Simulated, not measured. Final PBE-D3 vc-relax geometries (predicted structures); labels: hkl of the six strongest peaks.", ha="center", va="top", fontsize=6.5, color="#444")
    with open(HERE / "xrd_peaks.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["two_theta_deg", "d_A", "hkl", "rel_intensity", "structure"]); w.writerows(rows)
    save(fig, "F5_xrd", dy=-0.05)


# ---------------------------------------------------------------- T1 method-error / dispersion table
def fmt(x, d=2):
    return f"{x:+.{d}f}" if isinstance(x, float) else ("pending" if x is None else "not computed")


def t1():
    import csv
    tab = []   # (group, quantity, {method: value}, {method: [jobs]})

    def add(group, q, getter):
        vals, srcs = {}, {}
        for m in METHODS: vals[m], srcs[m] = getter(m)
        tab.append((group, q, vals, srcs))

    def margin(key, field):
        def g(m):
            v = tp_row(m, key)
            if v is None: return "not computed", []
            if v["status"] != "done": return None, [v["job"]]
            x = v[field]
            if not isinstance(x, float): return "not computed", [v["job"]]
            if field == "dE": return x, [v["job"], BEST[m][BASE][v["formula"]][2]]
            return x, [v["job"]] + [BEST[m][BASE][b][2] for b in v["binaries_route"]]
        return g

    def vs_libr(poly_source):
        def g(m):
            tpr, sb, lb = lowest(m, formula="LiSbBr4", tp=True), lowest(m, formula="SbBr3"), lowest(m, source=poly_source)
            if not (tpr and sb and lb): return ("pending" if any(r["method"] == m and r["status"] != "done" for r in ALL if r["formula"] == "LiBr" and r["source"] == poly_source) else "not computed"), []
            d = V.decompose("LiSbBr4", {"LiBr": lb["e_ry"] / lb["natoms"], "SbBr3": sb["e_ry"] / sb["natoms"]})
            return (tpr["e_ry"] / tpr["natoms"] - d[0]) * V.RY, [tpr["job"], lb["job"], sb["job"]]
        return g

    def order(src_exp, src_hyp):
        """E(experimental polymorph) - E(computationally lower competitor); positive = PBE-D3 ranks them in the wrong order."""
        def g(m):
            a, b = lowest(m, source=src_exp), lowest(m, source=src_hyp)
            if not (a and b): return "not computed", []
            return (a["e_ry"] / a["natoms"] - b["e_ry"] / b["natoms"]) * V.RY, [a["job"], b["job"]]
        return g

    add("A. Margin of the proposed model (negative = TP lower)", "TP-NaSbBr4-01 vs NaBr + lowest computed SbBr3", margin("Na", "binaries"))
    add("A", "TP-NaSbBr4-01 vs lowest listed NaSbBr4 (agm003733885)", margin("Na", "dE"))
    add("A", "TP-LiSbBr4-01 vs LiBr + SbBr3 (lowest LiBr computed in each column)", margin("Li", "binaries"))
    add("A", "TP-LiSbBr4-01 vs LiBr rocksalt (agm003157564) + SbBr3", vs_libr("agm003157564"))
    add("A", "TP-LiSbBr4-01 vs LiBr wurtzite (agm002091979) + SbBr3", vs_libr("agm002091979"))
    add("A", "TP-LiSbBr4-01 vs lowest listed LiSbBr4 (agm005739369)", margin("Li", "dE"))
    add("B. Polymorph-ordering error of the method (positive = wrong order)", "SbBr3: Pnma (expt, agm003218591) minus R-3 (agm002139075)", order("agm003218591", "agm002139075"))
    add("B", "LiBr: rocksalt (expt, agm003157564) minus wurtzite (agm002091979)", order("agm003157564", "agm002091979"))

    # reference check against numbers quoted in forensics_oct18 / chemrxiv.md (warn only)
    for (q_i, m, want) in ((6, "zero", 15.2), (6, "bj", 26.1), (6, "none", 20.6), (2, "zero", -7.66), (4, "none", -7.4)):
        got = tab[q_i][2][m] if q_i != 4 else tab[3][2][m]
        if not isinstance(got, float) or abs(got - want) > 0.15: print(f"NOTE: T1 row {q_i} {m} = {got} vs Oct-18 reference {want}")

    cols = ["group", "quantity (meV/atom)", *(f"{METHODS[m]}" for m in METHODS)]
    with open(DATA / "T1_method_error.csv", "w", newline="") as fh:
        fh.write(f"# {STAMP}; QE, {BASE}; meV/atom; 'pending' = job queued/unfinished, 'not computed' = no such job\n")
        w = csv.writer(fh); w.writerow(cols + [f"source jobs ({m})" for m in METHODS])
        for g, q, vals, srcs in tab:
            w.writerow([g[:1], q] + [fmt(vals[m], 6) for m in METHODS] + [" ; ".join(PATHS[j].removesuffix(".in") + ".out" for j in srcs[m]) for m in METHODS])
    L = [f"<!-- {STAMP}; generated by paper/figures/make_figures.py; do not edit -->", "",
         "# T1. Dispersion sensitivity of the proposed-model margins and PBE ordering errors of reference polymorphs", "",
         f"Quantum ESPRESSO, {BASE.replace('k', 'k-spacing ')}, 0 K, no SOC, meV/atom. Part A: negative = proposed model below the reference. "
         "Part B: energy of the experimentally observed polymorph minus the polymorph the method ranks lower; positive = wrong order. "
         "Predicted, proposed models; not an experimental claim.", ""]
    grp = None
    for g, q, vals, srcs in tab:
        if g != grp and len(g) > 1:
            if grp: L.append("")
            L += [f"**{g}**", "", "| Quantity | " + " | ".join(METHODS.values()) + " |", "|---|---|---|---|"]; grp = g
        L.append(f"| {q} | " + " | ".join(fmt(vals[m]) for m in METHODS) + " |")
    L += ["", "## Source job files (repo-relative; energy = last `!` line of the .out)", ""]
    for g, q, vals, srcs in tab:
        L.append(f"- **{q}**")
        for m in METHODS:
            L.append(f"  - {METHODS[m]}: {fmt(vals[m])}; " + (", ".join(f"`{PATHS[j].removesuffix('.in')}.out`" for j in srcs[m]) if srcs[m] else "no job"))
    (DATA / "T1_method_error.md").write_text("\n".join(L) + "\n")
    print("wrote T1_method_error.csv/.md")
    return tab


# ---------------------------------------------------------------- T2 structures + F4 + CIFs
CN_CUT = 3.6   # Angstrom: one fixed Br-neighbour cutoff for every coordination number (no CrystalNN); the models have no Sb/A-Br contact between 3.1 and 4.2 A
T2_ROWS = [("Na", "TP", dict(formula="NaSbBr4", tp=True)), ("Na", "listed", dict(source="agm003733885")),
           ("Li", "TP", dict(formula="LiSbBr4", tp=True)), ("Li", "listed1", dict(source="agm005739369")), ("Li", "listed2", dict(source="agm003740876"))]
SBBR3 = {"sbbr3_r3": ("R-3 SbBr$_3$ (agm002139075)", dict(source="agm002139075")),
         "sbbr3_pnma": ("Pnma SbBr$_3$ (agm003218591)", dict(source="agm003218591"))}


def _groups(d, gap=0.12):
    """Sorted distances -> '2+2+2' (a new group starts where consecutive distances differ by more than gap A)."""
    g = [[d[0]]]
    for x in d[1:]: (g[-1].append(x) if x - g[-1][-1] <= gap else g.append([x]))
    return "+".join(str(len(x)) for x in g)


def _dists(s, el, nn="Br", cut=CN_CUT):
    return [sorted(n.nn_distance for n in s.get_neighbors(site, cut) if n.specie.symbol == nn) for site in s if site.specie.symbol == el]


def _connectivity(s, cut=CN_CUT):
    """(Larsen dimensionality of the Sb-Br network, sorted set of shared-Br counts between Sb pairs: 1 corner, 2 edge, 3 face)."""
    from pymatgen.analysis.dimensionality import get_dimensionality_larsen
    from pymatgen.analysis.graphs import StructureGraph
    from pymatgen.analysis.local_env import CutOffDictNN
    dim = get_dimensionality_larsen(StructureGraph.from_local_env_strategy(s, CutOffDictNN({("Sb", "Br"): cut, ("Br", "Sb"): cut})))
    sb = [i for i, x in enumerate(s) if x.specie.symbol == "Sb"]
    nb = {i: {(n.index, tuple(int(v) for v in n.image)) for n in s.get_neighbors(s[i], cut) if n.specie.symbol == "Br"} for i in sb}
    pairs = set()
    for i in sb:
        for j in sb:
            for o in np.ndindex(3, 3, 3):
                o = np.array(o) - 1
                if i == j and not o.any(): continue
                c = len(nb[i] & {(k, tuple(int(v) for v in np.array(t) + o)) for k, t in nb[j]})
                if c: pairs.add(c)
    return dim, sorted(pairs)


def _conn_text(dim, pairs):
    share = {1: "corner", 2: "edge", 3: "face"}
    return f"{'molecular (0D)' if dim == 0 else str(dim) + 'D'} Sb-Br network, " + (" and ".join(share.get(p, f"{p}-Br") + "-sharing" for p in pairs) or "isolated") + " SbBr_n units"


def _rms_ang(a, b):
    """Loosest-tolerance ladder; RMS in Angstrom = pymatgen normalised rms * (V/N)^(1/3) of the volume-averaged primitive cell."""
    for tol in ((0.2, 0.3, 5), (0.3, 0.5, 10), (0.5, 0.7, 15), (0.7, 1.0, 20), (1.0, 1.5, 30)):
        r = StructureMatcher(*tol, primitive_cell=True, scale=True).get_rms_dist(a, b)
        if r is not None:
            pa, pb = a.get_primitive_structure(), b.get_primitive_structure(); f = ((pa.volume + pb.volume) / 2 / len(pa)) ** (1 / 3)
            return dict(rms_A=round(r[0] * f, 3), max_A=round(r[1] * f, 3), ltol_stol_angle_tol=list(tol))
    return dict(rms_A=None, max_A=None, ltol_stol_angle_tol="no match up to 1.0/1.5/30")


def t2_structures():
    import json
    cif = HERE / "cif"; cif.mkdir(exist_ok=True)
    S = {(k, role): (lowest("zero", **kw), None) for k, role, kw in T2_ROWS}
    S = {kr: (r, structure_of(r)) for kr, (r, _) in S.items()}
    rows, J = [], dict(cn_cutoff_A=CN_CUT, cn_method=f"fixed {CN_CUT} A cutoff on Br neighbours (not CrystalNN); PBE-D3 zero-damping relaxed cells, {BASE}", structures={})
    for (key, role), (r, s) in S.items():
        name = TPS[key]["tp"] if role == "TP" else r["source"]
        sga = SA(s, 0.05); ref = sga.get_refined_structure(); sym = SA(ref, 0.05).get_symmetrized_structure(); L = ref.lattice
        z = round(ref.composition.get_reduced_composition_and_factor()[1]); sg, n = sga.get_space_group_symbol(), sga.get_space_group_number()
        wy = "; ".join(f"{w} {e[0].specie.symbol}: " + " ".join(f"{x:.3f}" for x in e[0].frac_coords) for w, e in zip(sym.wyckoff_symbols, sym.equivalent_sites))
        rows.append(f"| {name} ({r['formula']}) | {sg} ({n}) | {L.a:.3f} | {L.b:.3f} | {L.c:.3f} | {L.alpha:.1f} / {L.beta:.2f} / {L.gamma:.1f} | {s.volume / (len(s) / 6):.1f} | {z} | {wy} |")
        CifWriter(s, symprec=0.05).write_file(cif / f"{name}_{r['formula']}_PBE-D3_relaxed.cif")
        J["structures"][name] = dict(job=r["job"], space_group=sg, number=n, a=L.a, b=L.b, c=L.c, alpha=L.alpha, beta=L.beta, gamma=L.gamma, V_per_fu_A3=s.volume / (len(s) / 6), Z=z)
    for key in TPS:
        s = S[(key, "TP")][1]; d = _dists(s, "Sb")[0]; a = _dists(s, TPS[key]["f"][:2])[0]; dim, pairs = _connectivity(s); cn = len(d)
        J[key.lower()] = dict(sb_cn=cn, sbbr_min=round(min(d), 3), sbbr_max=round(max(d), 3), a_cn=len(a), a_br_range_A=[round(min(a), 3), round(max(a), 3)],
                              sb_motif=f"{'distorted octahedral' if cn == 6 else str(cn) + '-fold'} ({_groups(d)})", sb_connectivity=_conn_text(dim, pairs),
                              sb_dimensionality=dim, shared_br_counts=pairs, sbbr_all_A=[round(x, 3) for x in d])
    for k, (lab, kw) in SBBR3.items():
        s = structure_of(lowest("zero", **kw)); d = _dists(s, "Sb")[0]; dim, pairs = _connectivity(s); cn = len(d)
        outer = [round(x, 2) for x in _dists(s, "Sb", cut=4.2)[0] if x > CN_CUT]
        J[k] = dict(sb_cn=cn, sbbr_min=round(min(d), 3), sbbr_max=round(max(d), 3), sb_dimensionality=dim, shared_br_counts=pairs, secondary_contacts_3p6_to_4p2_A=outer,
                    sb_motif=f"{ {3: 'trigonal-pyramidal', 6: 'distorted octahedral'}.get(cn, str(cn) + '-fold') } ({_groups(d)})",
                    sb_connectivity=(f"isolated SbBr_3 molecules (secondary Sb...Br contacts at {outer[0]}-{outer[-1]} A, beyond the {CN_CUT} A cutoff)" if dim == 0 and outer else _conn_text(dim, pairs)))
    r3, na, li = J["sbbr3_r3"], J["na"], J["li"]
    share = na["sb_cn"] == li["sb_cn"] == r3["sb_cn"] == 6 and min(na["sb_dimensionality"], li["sb_dimensionality"], r3["sb_dimensionality"]) > 0
    J["shared_motif_vs_r3"] = (f"six-coordinate, distorted-octahedral Sb (R-3: {r3['sb_motif']}; models: {na['sb_motif']} and {li['sb_motif']}) in an extended Sb-Br network"
                               if share else "no common motif")
    J["na_vs_listed"] = _rms_ang(S[("Na", "TP")][1], S[("Na", "listed")][1]) | dict(listed="agm003733885")
    J["li_vs_listed2"] = _rms_ang(S[("Li", "TP")][1], S[("Li", "listed2")][1]) | dict(listed="agm003740876")
    J["li_vs_listed1"] = _rms_ang(S[("Li", "TP")][1], S[("Li", "listed1")][1]) | dict(listed="agm005739369")
    a, b = structure_of(lowest("zero", source="agm003218590")), structure_of(lowest("zero", source="agm003218591"))
    J["sbbr3_p212121_start"] = dict(symprec_0p001=SA(a, 0.001).get_space_group_symbol(), symprec_0p05=SA(a, 0.05).get_space_group_symbol(),
                                    pnma_start_symprec_0p001=SA(b, 0.001).get_space_group_symbol(), volume_A3=[round(a.volume, 1), round(b.volume, 1)],
                                    structurematcher_default_match=bool(StructureMatcher().fit(a, b)),
                                    note="relaxed P2_12_12_1-start cell stays P2_12_12_1 at symprec 0.001 and only reads as Pnma at 0.05; not the same structure as the Pnma-start cell")
    (DATA / "T2_structures.json").write_text(json.dumps(J, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    hdr = ["| Structure | Space group (No.) | a (Å) | b (Å) | c (Å) | α / β / γ (°) | V/f.u. (Å^3) | Z | Wyckoff positions (site: x y z) |", "|---|---|---|---|---|---|---|---|---|"]
    (DATA / "T2_structures.md").write_text(f"# {STAMP}; conventional (refined) cells, symprec 0.05; relaxed PBE-D3 zero damping, {BASE}\n" + "\n".join(hdr + rows) + "\n")
    print("wrote T2_structures.md/.json and", len(rows), "CIFs"); return S, J


def _draw_cell(ax, s):
    """Sb-Br polyhedra (convex hull of bonded Br) and A cations in a 2x2x2 supercell; depth-sorted mplot3d faces."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from scipy.spatial import ConvexHull
    s = SpacegroupAnalyzer(s, symprec=0.05).get_conventional_standard_structure()  # axes as in Table 2
    sc = s.copy(); sc.make_supercell((2, 2, 2)); ctr = sc.lattice.get_cartesian_coords([.5, .5, .5]); R = 0.55 * np.linalg.norm(s.lattice.matrix.sum(0))
    near = lambda p: np.linalg.norm(p - ctr) < R
    brs = []
    for site in sc:
        el = site.specie.symbol
        if not near(site.coords): continue
        if el == "Sb":
            nb = np.array([n.coords for n in sc.get_neighbors(site, CN_CUT) if n.specie.symbol == "Br"])
            ax.add_collection3d(Poly3DCollection([nb[f] for f in ConvexHull(nb).simplices], facecolor=C_TP, alpha=0.28, edgecolor="#7a3500", lw=0.4))
            ax.scatter(*site.coords, s=24, color=C_TP, edgecolor="k", lw=0.3, depthshade=False); brs += [tuple(np.round(x, 3)) for x in nb]
        elif el in ("Na", "Li"): ax.scatter(*site.coords, s=26, color=C_LISTED, edgecolor="k", lw=0.3, depthshade=False)
    if brs: ax.scatter(*np.array(sorted(set(brs))).T, s=12, color="#8C6D31", depthshade=False)  # shared Br drawn once
    o = ctr - s.lattice.matrix.sum(0) / 2
    for v, nm in zip(s.lattice.matrix, "abc"): ax.plot(*zip(o, o + v), color="k", lw=0.8); ax.text(*(o + v * 1.1), nm, fontsize=7)
    r = 0.72 * R; ax.set_xlim(ctr[0] - r, ctr[0] + r); ax.set_ylim(ctr[1] - r, ctr[1] + r); ax.set_zlim(ctr[2] - r, ctr[2] + r)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(elev=20, azim=-60)


def f4(S, J):
    fig = plt.figure(figsize=(7.2, 8.2)); gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.1], hspace=0.10)
    top = gs[0].subgridspec(1, 2, wspace=0.0); bot = gs[1].subgridspec(4, 1, hspace=0.15)
    for j, key in enumerate(TPS):
        ax = fig.add_subplot(top[j], projection="3d"); _draw_cell(ax, S[(key, "TP")][1]); tp = TPS[key]["tp"]
        ax.set_title(f"({'ab'[j]}) {tp}, {J['structures'][tp]['space_group']}\nSb: {J[key.lower()]['sb_motif']}", loc="left", fontsize=7, fontweight="bold")
    sets = [(TPS["Na"]["tp"], S[("Na", "TP")][1], C_TP), (TPS["Li"]["tp"], S[("Li", "TP")][1], C_TP)] + \
           [(lab, structure_of(lowest("zero", **kw)), C_LISTED) for lab, kw in SBBR3.values()]
    bins = np.arange(2.4, 4.05, 0.05); axs = []
    for k, (lab, s, c) in enumerate(sets):
        ax = fig.add_subplot(bot[k], sharex=axs[0] if axs else None); axs.append(ax); dd = _dists(s, "Sb", cut=4.0)
        ax.hist([x for d in dd for x in d], bins=bins, weights=np.full(sum(map(len, dd)), 1 / len(dd)), color=c, edgecolor="k", lw=0.3)
        ax.axvline(CN_CUT, color=C_GREY, ls="--", lw=0.7); ax.set_ylabel("Br per Sb", fontsize=6.5)
        ax.text(0.2 if k == 3 else 0.45, 0.92, lab, transform=ax.transAxes, va="top", fontsize=7, fontweight="bold")
        if k < 3: plt.setp(ax.get_xticklabels(), visible=False)
    axs[0].text(CN_CUT + 0.02, 0.5, f"CN cutoff {CN_CUT} Å", fontsize=6, color="#555", transform=axs[0].get_xaxis_transform())
    axs[-1].set_xlabel("Sb–Br distance (Å)")
    fig.text(0.5, 0.02, "(a, b) orange: SbBr$_n$ polyhedra (Sb–Br up to 3.6 Å); brown: Br; blue: Na/Li; conventional C2/c cells, chains along c. (c) Sb–Br distance histograms, 0.05 Å bins. PBE-D3 (zero damping) relaxed cells. Proposed models, predicted.",
             ha="center", fontsize=6.3, color="#444")
    save(fig, "F4_structures", dy=0.0)


if __name__ == "__main__":
    print(STAMP)
    f1(); f2(); f3(); f5(); t1(); S2, J2 = t2_structures(); f4(S2, J2)
