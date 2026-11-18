"""M5 nights 6+ (from 2026-10-14). Run from repo root:  discovery/.venv/bin/python -I overnight/make_nights6plus.py [nightN]
1. Writes the vc-relax inputs in PLAN (missing references / competitors; settings via make_inputs.write: PBE-D3, PSlibrary PAW,
   60/480 Ry where compared with the NaBiCl4/NaSbCl4 60 Ry results, 80/640 Ry for Li and bromide systems, k-spacing 0.22).
2. G4 step: for every FINISHED vc-relax in a candidate's chemical system at the candidate's base setting that has no
   80 Ry / k-spacing 0.17 SCF anywhere under overnight/*/jobs, writes one (make_converge_nasbbr4.py approach: same input,
   final geometry, calculation='scf', 80/640 Ry, denser k) into nightN (default night7). Re-run after every morning push
   with the next unrun night as argument; it never touches an input that already exists and appends new nights to QUEUE.
Starting geometries: MACE+D3 cells already in discovery/results/*/relax when present, else the Alexandria PBE cell
(cached in overnight/start/<id>.cif). Experimental (MP-flagged) structures are referenced by the Alexandria ID that
StructureMatcher maps them to (default tolerances, 2026-10-14). No structure is invented; every ID is an Alexandria entry."""
import glob, importlib.util, json, re, sys, urllib.parse, urllib.request
from pathlib import Path
from pymatgen.core import Composition, Structure
O = Path(__file__).resolve().parent; ROOT = O.parent
sys.path.insert(0, str(O))
import make_inputs as mi
import run_night as rn
_s = importlib.util.spec_from_file_location("verdict", ROOT / "discovery/scripts/verdict.py"); V = importlib.util.module_from_spec(_s); _s.loader.exec_module(V)

# night -> [(job name without number, Alexandria ID, ecut)]; est. M5 wall in the comment; what it serves after "->"
PLAN = {
 "night6": [("LiCl_rocksalt_agm003157566", 80),      # 2 at, 0.1 h -> G3 LiInCl4, LiSbCl4: experimental LiCl (Fm-3m) never computed; "LiCl_rs" jobs are zincblende
            ("LiBr_rocksalt_agm003157564", 80),      # 2 at, 0.1 h -> G3 LiSbBr4, LiBiBr4, LiScBr4: same gap for LiBr
            ("Na3SbCl6_agm067283895", 60),           # 10 at, 0.7 h -> G3 NaSbCl4: MACE+D3 9 meV/atom below the computed agm017656109
            ("Na3SbCl6_MPexpt_agm005507662", 60),    # 10 at, 0.5 h -> G3 NaSbCl4 caveat: = Materials Project Na3SbCl6 (Fm-3m)
            ("SbCl3_Pnma_expt_agm003219609", 60),    # 16 at, 1.2 h -> G3 NaSbCl4: experimental SbCl3 (MP Pnma)
            ("BiCl3_agm002142653", 60),              # 8 at, 0.5 h -> C11 NaBiCl4: second BiCl3 polymorph (Alexandria hull 0)
            ("BiCl3_Pnma_expt_agm003218787", 60),    # 16 at, 1.2 h -> C11 NaBiCl4: experimental BiCl3 (MP Pnma)
            ("SbCl3_agm002205242", 80),              # 8 at, 0.5 h -> G3 LiSbCl4: lowest SbCl3 at 60 Ry, no 80 Ry relax yet
            ("SbBr3_agm002257875", 80),              # 8 at, 0.5 h -> G3 LiSbBr4, NaSbBr4: second SbBr3 polymorph
            ("SbBr3_Pnma_expt_agm003218591", 80),    # 16 at, 1.3 h -> G3 LiSbBr4, NaSbBr4: experimental SbBr3 (MP Pnma, MP hull 0)
            ("Na3SbCl6_agm003569162", 60),           # 20 at, 2.0 h -> G3 NaSbCl4: MACE+D3 also ~9 meV/atom below agm017656109
            ("BiBr3_C2m_expt_agm002330398", 80)],    # 8 at, 0.6 h -> G3 LiBiBr4: experimental BiBr3 (MP C2/m)
 "night7": [("ScBr3_agm002184236", 80),              # 8 at, 0.6 h -> G3 LiScBr4: R-3 polymorph (Alexandria hull 0; C2/m one is on the Air)
            ("InCl3_agm002142775", 80),              # 8 at, 0.6 h -> G3 LiInCl4: R-3 polymorph (Alexandria hull 0)
            ("SbCl3_Pnma_expt_agm003219609", 80),    # 16 at, 1.3 h -> G3 LiSbCl4: experimental SbCl3 at 80 Ry
            ("SbBr3_P212121_expt_agm003218590", 80)],# 16 at, 1.3 h -> G3 LiSbBr4, NaSbBr4: alpha-SbBr3 (MP P2_12_12_1)
 "night8": [("Li2SbCl5_agm005868832_ternary", 80),   # 16 at, 1.3 h -> G3 LiSbCl4: Li-Sb-Cl ternary nearest the hull (+41 PBE)
            ("Na2SbCl5_agm005868900_ternary", 60),   # 8 at, 0.5 h -> G3 NaSbCl4: Na-Sb-Cl ternary on the NaCl-SbCl3 line (+26 PBE)
            ("Na2SbBr5_agm005863603_ternary", 80),   # 8 at, 0.5 h -> G3 NaSbBr4: same for Na-Sb-Br (+31 PBE)
            ("Li3SbBr6_P4332_agm070912415_ternary", 80)],  # 40 at, 5-8 h -> G3 LiSbBr4: lowest listed Li3SbBr6 (45 meV/atom below the computed 10-atom cell in PBE)
 "night9": [("Li3BiBr6_P4332_agm071391430_ternary", 80)],  # 40 at, 5-8 h -> G3 LiBiBr4: same for Li3BiBr6
}
# candidate chemical systems still alive (verdict.py 2026-10-14; TP-LiBiCl4-01 failed G3) -> base setting
SYSTEMS = {"Bi Cl Na": "60 Ry k0.22", "Cl Na Sb": "60 Ry k0.22", "Br Na Sb": "80 Ry k0.22", "Br Li Sb": "80 Ry k0.22",
           "Bi Br Li": "80 Ry k0.22", "Cl Li Sb": "80 Ry k0.22", "Cl In Li": "80 Ry k0.22", "Br Li Sc": "80 Ry k0.22",
           "Bi Br Na": "80 Ry k0.22", "I Na Sb": "80 Ry k0.22",  # 2026-10-15 new candidates (new_candidates_oct15.md)
           "Br Li Y": "80 Ry k0.22", "Ga I Na": "80 Ry k0.22"}  # 2026-10-16: LiYBr4 passes G3, NaGaI4 G2 pass


def start(aid):
    for d in ("gate5", "gate4", "gate3", "skeptic3"):
        p = ROOT / f"discovery/results/{d}/relax/{aid}_d3.cif"
        if p.exists(): return Structure.from_file(p)
    p = O / "start" / f"{aid}.cif"
    if not p.exists():
        u = "https://alexandria.icams.rub.de/pbe/v1/structures?" + urllib.parse.urlencode({"filter": f'id="{aid}"',
            "response_fields": "lattice_vectors,species_at_sites,cartesian_site_positions"})
        a = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "halide/0.1"}), timeout=120))["data"][0]["attributes"]
        p.parent.mkdir(exist_ok=True)
        Structure(a["lattice_vectors"], a["species_at_sites"], a["cartesian_site_positions"], coords_are_cartesian=True).to(filename=str(p))
    return Structure.from_file(p)


def queue(night):
    q = O / "QUEUE"; names = q.read_text().split()
    if night not in names: q.write_text("\n".join(names + [night]) + "\n")


def nxt(jobdir):
    return max([int(m.group(1)) for f in jobdir.glob("*.in") if (m := re.match(r"(\d+)_", f.name))] or [0]) + 1


def relaxes():
    for night, jobs in PLAN.items():
        mi.OUT = O / night / "jobs"; mi.OUT.mkdir(parents=True, exist_ok=True); queue(night)
        for n, (name, ecut) in enumerate(jobs, 1):
            full = f"{n:02d}_{name}_vcrelax"
            if not (mi.OUT / f"{full}.in").exists():
                mi.write(full, start(re.search(r"agm\d+", name).group(0)), "vc-relax", ecut=float(ecut))


def converge(night):
    jobs = [Path(f) for f in sorted(glob.glob(str(O / "*/jobs/*.in"))) if Path(f).parent.parent.name not in V.SKIP]
    have = []  # (formula, natoms, volume/atom) of every 80 Ry k0.17 SCF input, queued or done (same geometry => same volume)
    for f in jobs:
        if "_k017" in f.stem and "'scf'" in f.read_text():
            s = mi.final_structure(f); have.append((s.composition.reduced_formula, len(s), s.volume / len(s)))
    out_dir = O / night / "jobs"; out_dir.mkdir(parents=True, exist_ok=True); n = nxt(out_dir); new = []
    for f in jobs:
        r = V.parse(f)
        if r["status"] != "done" or r["calc"] != "vc-relax": continue
        els = {e.symbol for e in Composition(r["formula"])}
        if not any(els <= set(k.split()) and r["setting"] == st for k, st in SYSTEMS.items()): continue
        out = f.with_suffix(".out"); s = mi.final_structure(out)
        if any(fm == r["formula"] and na == len(s) and abs(v - s.volume / len(s)) < 1e-3 for fm, na, v in have): continue
        stem = re.sub(r"_(60|80)Ry", "", re.sub(r"^\w+?_", "", f.stem, count=1).removesuffix("_vcrelax"))
        job = f"{n:02d}_{stem}_80Ry_k017_scf"; n += 1
        txt = rn.with_geometry(f.read_text(), rn.last_geometry(out.read_text(errors="replace")))
        txt = txt.replace("calculation='vc-relax'", "calculation='scf'").replace("ecutwfc=60.0", "ecutwfc=80.0").replace("ecutrho=480.0", "ecutrho=640.0")
        txt = re.sub(r"prefix='[^']*'", f"prefix='{job}'", txt)
        txt = re.sub(r"outdir='__TMP__/[^']*'", f"outdir='__TMP__/{job}'", txt)
        k = mi.kmesh(s, 0.17)
        txt = re.sub(r"(K_POINTS automatic\n)\s*\d+ \d+ \d+", lambda m: m.group(1) + "  " + " ".join(map(str, k)), txt)
        assert "ecutwfc=80.0" in txt and "'scf'" in txt, f
        (out_dir / f"{job}.in").write_text(txt); have.append((r["formula"], len(s), s.volume / len(s))); new.append(job)
        print(job, "from", r["job"], "k", k)
    if new: queue(night)
    return new


if __name__ == "__main__":
    # The PLAN relaxations were written on 2026-10-14 and then moved into date-named batches (overnight/m5_oct14..18);
    # relaxes() is not re-run so it cannot recreate them under the old night6-9 names.
    assert len(sys.argv) == 2, "usage: make_nights6plus.py <batch, e.g. m5_oct16>  (queues G4 SCFs for newly finished relaxations)"
    converge(sys.argv[1])
