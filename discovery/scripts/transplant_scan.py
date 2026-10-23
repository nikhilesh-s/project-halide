"""Cycle 3: prototype-transplant scan for A M X4 halides (A=Li,Na). Everything produced is a PROPOSED MODEL
ranked by MACE-MPA-0 (+D3); nothing here is a stability or discovery claim.
Resumable (one json per relaxation, downloads cached). Run from repo root:
    discovery/.venv/bin/python -I discovery/scripts/transplant_scan.py [cheap_minutes=25] [d3_minutes=20] [max_proto_atoms=24] [need_nonF=0]
Scope shrinks by priority order if a stage deadline hits: iodides dropped first, then Sc/Y.
max_proto_atoms / need_nonF (prototype group must contain a Cl/Br/I member) cut the prototype list further; rerun with
bigger values to extend the scan (finished relaxations are reused)."""
import json, sys, time, csv, shutil, urllib.request, urllib.parse
from pathlib import Path
from multiprocessing import Pool
import numpy as np
from pymatgen.core import Structure, Lattice, Element
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

D = Path("discovery"); CACHE = D/"cache"; RUN = D/"results"/"transplant_cycle3"; CAND = D/"candidates"
U = "https://alexandria.icams.rub.de/pbe/v1/structures"
A_T, M_ALL, X_T = ["Na", "Li"], ["Bi", "Sb", "Al", "Ga", "In", "Sc", "Y"], ["Cl", "Br", "I"]
A_P, X_P = ["Li", "Na", "K", "Rb", "Cs", "Ag", "Tl"], ["F", "Cl", "Br", "I"]
POOL_HULL, MAXN, KEEP_EV, HIT_EV = 0.030, 24, 0.015, 0.003
NPROC = 8
SM = StructureMatcher()                                    # 0.2 / 0.3 / 5 deg
SML = StructureMatcher(ltol=0.3, stol=0.5, angle_tol=10)   # "loose", as in skeptic/identity.py

def rform(d):  # OPTIMADE reduced formula: alphabetical, count omitted when 1
    return "".join(e + (str(n) if n > 1 else "") for e, n in sorted(d.items()))

def fetch(form):
    f = CACHE/f"alex_{form}.json"
    if f.exists(): return json.load(open(f))
    url = U + "?" + urllib.parse.urlencode({"filter": f'chemical_formula_reduced="{form}"', "page_limit": 100,
          "response_fields": "lattice_vectors,cartesian_site_positions,species_at_sites,nsites,_alexandria_hull_distance,_alexandria_space_group"})
    out = []
    while url:
        for k in range(5):
            try: r = json.load(urllib.request.urlopen(url, timeout=60)); break
            except Exception as e: err = e; time.sleep(3*(k+1))
        else: raise err
        out += [{"id": d["id"], **d["attributes"]} for d in r["data"]]
        nx = (r.get("links") or {}).get("next"); url = nx.get("href") if isinstance(nx, dict) else nx
        time.sleep(0.5)                                    # <= 2 requests / s
    json.dump(out, open(f, "w")); return out

def struct(e):
    return Structure(Lattice(e["lattice_vectors"]), e["species_at_sites"], e["cartesian_site_positions"], coords_are_cartesian=True)

# ---------- relaxation workers (settings identical to skeptic/run.py; cheap pass is looser by design) ----------
_C = None
def _init(d3):   # note: cheap pass runs float32 (1.6x faster, ~0.1 meV/atom noise); D3 pass is float64 as in the skeptic scripts
    global _C
    import torch; torch.set_num_threads(1)
    from mace.calculators import mace_mp
    _C = mace_mp(model="medium-mpa-0", device="cpu", default_dtype="float64" if d3 else "float32", dispersion=d3)   # D3: BJ damping, PBE params (mace_mp default)

def _relax(job):
    name, sd, fmax, steps, out, deadline = job
    f = Path(out)/f"{name}.json"
    if f.exists(): return json.load(open(f))
    if time.time() > deadline: return None
    from pymatgen.io.ase import AseAtomsAdaptor as AA
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE
    s0 = Structure.from_dict(sd); a = AA.get_atoms(s0); a.calc = _C; n = len(a); t = time.time()
    try:
        opt = FIRE(FrechetCellFilter(a), logfile=None)
        conv = bool(opt.run(fmax=fmax, steps=steps))
        s1 = AA.get_structure(a); dm = s1.distance_matrix + np.eye(n)*1e6
        r = dict(name=name, n=n, E=float(a.get_potential_energy()/n), V=s1.volume/n, V0=s0.volume/n, conv=conv,
                 steps=int(opt.nsteps), min_dist=float(dm.min()), wall=time.time()-t, structure=s1.as_dict())
    except Exception as e:
        r = dict(name=name, n=n, error=repr(e)[:200], wall=time.time()-t)
    json.dump(r, open(f, "w")); return r

def run_stage(jobs, d3, tag):
    out = RUN/tag; out.mkdir(parents=True, exist_ok=True)
    todo = [j for j in jobs if not (out/f"{j[0]}.json").exists()]
    print(f"[{tag}] {len(jobs)} jobs, {len(todo)} to run", flush=True)
    if todo:
        with Pool(NPROC, initializer=_init, initargs=(d3,)) as p:
            for i, r in enumerate(p.imap_unordered(_relax, todo, chunksize=1)):
                if i % 100 == 0: print(f"  [{tag}] {i}/{len(todo)} t={time.time()-T0:.0f}s", flush=True)
    return {j[0]: json.load(open(out/f"{j[0]}.json")) for j in jobs if (out/f"{j[0]}.json").exists()}

def ok(r): return r is not None and "E" in r and r["min_dist"] > 1.5
def prim(s): return s.get_primitive_structure()
def same(a, b):
    """'default' / 'loose' / '' : species-aware match of two primitive structures."""
    if len(a) != len(b): return ""
    return "default" if SM.fit(a, b) else "loose" if SML.fit(a, b) else ""
def sgs(s):
    try: return SpacegroupAnalyzer(s, symprec=0.01).get_space_group_symbol()
    except Exception: return "?"

if __name__ == "__main__":
    T0 = time.time()
    cheap_min = float(sys.argv[1]) if len(sys.argv) > 1 else 25
    d3_min = float(sys.argv[2]) if len(sys.argv) > 2 else 20
    MAXP = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    NONF = len(sys.argv) > 4 and sys.argv[4] == "1"
    for p in (CACHE, RUN): p.mkdir(parents=True, exist_ok=True)

    # ---- 1. downloads: all 196 A'M'X'4 chemistries (targets are a subset) ----
    ENT = {}
    for A in A_P:
        for M in M_ALL:
            for X in X_P:
                es = [e for e in fetch(rform({A: 1, M: 1, X: 4})) if e["nsites"] <= MAXN]
                for e in es: e.update(A=A, M=M, X=X)
                ENT[A, M, X] = es
    print("fetched", sum(map(len, ENT.values())), "entries <=24 atoms in", len(ENT), "chemistries", f"t={time.time()-T0:.0f}s", flush=True)
    BYID = {e["id"]: e for es in ENT.values() for e in es}
    def dummy(e):  # role-preserving decoration A->Na, M->Bi, X->Cl so species-aware matching compares roles
        s = struct(e); s.replace_species({e["A"]: "Na", e["M"]: "Bi", e["X"]: "Cl"}); return prim(s)
    DUM = {i: dummy(e) for i, e in BYID.items()}

    # ---- pool -> distinct prototypes. note: role-preserving fit (A/M/X roles fixed) instead of fit_anonymous;
    # it is the stricter test (an A<->M swapped twin counts as a separate prototype, which is what a transplant needs).
    pool = [e for e in BYID.values() if e["_alexandria_hull_distance"] <= POOL_HULL]
    gf = CACHE/"prototype_groups.json"
    if gf.exists(): GROUPS = json.load(open(gf))
    else:
        inv = {id(DUM[e["id"]]): e["id"] for e in pool}
        GROUPS = [sorted((inv[id(s)] for s in g), key=lambda i: BYID[i]["_alexandria_hull_distance"])
                  for g in SM.group_structures([DUM[e["id"]] for e in pool])]
        json.dump(GROUPS, open(gf, "w"))
    REPS = [g[0] for g in GROUPS]; GOF = {i: g[0] for g in GROUPS for i in g}
    print(f"pool {len(pool)} entries -> {len(REPS)} prototypes t={time.time()-T0:.0f}s", flush=True)

    # ---- 2. transplants ----
    prio = lambda c: (c[2] == "I", c[1] in ("Sc", "Y"), (c[0], c[1]) != ("Na", "Bi"))
    TARGETS = sorted(((A, M, X) for A in A_T for M in M_ALL for X in X_T), key=prio)
    def transplant(rep, c):
        e = BYID[rep]; s = struct(e); mp = {e["A"]: c[0], e["M"]: c[1], e["X"]: c[2]}
        v0 = s.volume/len(s); ex = ENT[c]
        if ex: v = float(np.mean([struct(x).volume/x["nsites"] for x in ex])); how = "mean V/atom of existing Alexandria entries"
        else:
            r3 = lambda a, m, x: Element(a).average_ionic_radius**3 + Element(m).average_ionic_radius**3 + 4*Element(x).average_ionic_radius**3
            v = v0*float(r3(*c)/r3(e["A"], e["M"], e["X"])); how = "sum of ionic-radius^3 ratio (no existing entries)"
        s.replace_species({k: w for k, w in mp.items() if k != w}); s.scale_lattice(v*len(s))
        return s, dict(source_id=rep, source_formula=f'{e["A"]}{e["M"]}{e["X"]}4', source_hull_eV=e["_alexandria_hull_distance"],
                       source_sg=e["_alexandria_space_group"], substitution=mp, volume_scale=v/v0, target_V_per_atom=v, scale_rule=how,
                       prototype_group=[i for i in GROUPS[REPS.index(rep)]])
    USE = {g[0] for g in GROUPS if BYID[g[0]]["nsites"] <= MAXP and (not NONF or any(BYID[i]["X"] != "F" for i in g))}
    print(f"prototypes in scope this run: {len(USE)} of {len(REPS)}", flush=True)
    ROWS, cheap_jobs, DL = {}, [], time.time() + cheap_min*60
    xf = CACHE/"excluded.json"; XC = json.load(open(xf)) if xf.exists() else {}
    for c in TARGETS:
        cn = "".join(c)
        if cn not in XC:
            XC[cn] = {}
            for rep in REPS:
                m = next((x["id"] for x in ENT[c] if GOF.get(x["id"]) == rep), None) or \
                    next((x["id"] for x in ENT[c] if len(DUM[x["id"]]) == len(DUM[rep]) and SM.fit(DUM[rep], DUM[x["id"]])), None)
                if m: XC[cn][rep] = m
        for x in ENT[c]: cheap_jobs.append((f"ex_{x['id']}", struct(x).as_dict(), 0.05, 150, str(RUN/"cheap"), DL))
        for rep in REPS:
            row = ROWS[cn, rep] = dict(chem=cn + "4", proto=rep, proto_formula="".join(BYID[rep][k] for k in "AMX") + "4",
                                       n_atoms=BYID[rep]["nsites"], n_existing=len(ENT[c]), status="")
            if rep in XC[cn]: row["status"] = "already_listed"; row["matches"] = XC[cn][rep]; continue
            if rep not in USE and not (RUN/"cheap"/f"tp_{cn}_{rep}.json").exists(): row["status"] = "skipped_prototype_scope"; continue
            s, prov = transplant(rep, c); row["volume_scale"] = round(prov["volume_scale"], 4)
            cheap_jobs.append((f"tp_{cn}_{rep}", s.as_dict(), 0.05, 150, str(RUN/"cheap"), DL))
    json.dump(XC, open(xf, "w"))
    n_att = sum(j[0].startswith("tp_") for j in cheap_jobs)
    print(f"transplants to attempt: {n_att}; existing entries: {len(cheap_jobs)-n_att} t={time.time()-T0:.0f}s", flush=True)
    CH = run_stage(cheap_jobs, False, "cheap")

    # ---- cheap-pass triage ----
    S = lambda r: prim(Structure.from_dict(r["structure"]))
    d3_jobs, DL, INFO = [], time.time() + d3_min*60, {}
    ANCH = {"NaBiBr": ["agm069058740"], "NaBiCl": ["agm088474086", "agm005743614"]}
    for c in TARGETS:
        cn = "".join(c); exn = [f"ex_{x['id']}" for x in ENT[c]]; tpn = [k for k in ROWS if k[0] == cn and ROWS[k]["status"] == ""]
        if any(n not in CH for n in exn) or any(f"tp_{cn}_{r}" not in CH for _, r in tpn):
            for k in tpn: ROWS[k]["status"] = "dropped_budget"
            continue
        exok = sorted((n for n in exn if ok(CH[n])), key=lambda n: CH[n]["E"])
        if not exok:
            for k in tpn: ROWS[k]["status"] = "no_existing_reference"
            continue
        best = CH[exok[0]]["E"]; exS = {n: [prim(struct(BYID[n[3:]]))] + ([S(CH[n])] if ok(CH[n]) else []) for n in exn}; kept = []
        for k in tpn:
            row = ROWS[k]; r = CH[f"tp_{cn}_{k[1]}"]
            if not ok(r): row["status"] = "cheap_failed"; continue
            row.update(dE_cheap_meV=round(1000*(r["E"]-best), 2), cheap_conv=r["conv"], V_cheap=round(r["V"], 2))
            if r["E"]-best > KEEP_EV: row["status"] = "cheap_above_window"; continue
            s = S(r); m = next(((n[3:], t) for n, ss in exS.items() for x in ss if (t := same(s, x))), None)
            if m: row.update(status="rediscovery_cheap", matches=f"{m[0]}({m[1]})"); continue
            d = next((kk for kk, ss in kept if abs(CH[f"tp_{cn}_{kk}"]["E"]-r["E"]) < 0.002 and same(s, ss) == "default"), None)
            if d: row.update(status="duplicate_of_other_transplant", matches=d); continue
            kept.append((k[1], s)); row["status"] = "kept_for_D3"
        ex3 = exok[:3] + [f"ex_{i}" for i in ANCH.get(cn, []) if f"ex_{i}" in exok and f"ex_{i}" not in exok[:3]]
        INFO[c] = dict(ex=exok, ex3=ex3, kept=[k for k, _ in kept], best_cheap=exok[0][3:])
        if kept or cn in ANCH:                       # anchor chemistries always get their D3 reference relaxations
            for n in ex3 + [f"tp_{cn}_{k}" for k, _ in kept]:
                d3_jobs.append((n, CH[n]["structure"], 0.01, 600, str(RUN/"d3"), DL))
    n_kept = sum(len(v["kept"]) for v in INFO.values())
    print(f"cheap pass complete for {len(INFO)}/{len(TARGETS)} chemistries; kept {n_kept} transplants t={time.time()-T0:.0f}s", flush=True)
    D3 = run_stage(d3_jobs, True, "d3")

    # ---- D3 ranking + identity re-check; then verify hit chemistries against ALL their existing entries with D3 ----
    def rank(c, exd3):
        cn = "".join(c); info = INFO[c]; exd = [n for n in exd3 if ok(D3.get(n))]
        if not exd: return []
        bn = min(exd, key=lambda n: D3[n]["E"]); best = D3[bn]["E"]; hits = []
        exS = {n: [prim(struct(BYID[n[3:]])), S(CH[n])] + ([S(D3[n])] if ok(D3.get(n)) else []) for n in info["ex"]}
        exS.update({f"ex_{x['id']}": [prim(struct(x))] for x in ENT[c] if f"ex_{x['id']}" not in exS})
        for rep in info["kept"]:
            row = ROWS[cn, rep]; r = D3.get(f"tp_{cn}_{rep}")
            if r is None: row["status"] = "dropped_budget_D3"; continue
            if not ok(r): row["status"] = "D3_failed"; continue
            s = S(r)
            row.update(dE_D3_meV=round(1000*(r["E"]-best), 2), best_existing_D3=bn[3:], D3_conv=r["conv"], V_D3=round(r["V"], 2), sg_D3=sgs(s),
                       n_existing_D3=len(exd))
            m = next(((n[3:], t) for n, ss in exS.items() for x in ss if (t := same(s, x))), None)
            if m: row.update(status="rediscovery_D3", matches=f"{m[0]}({m[1]})"); continue
            row["status"] = "survives_identity"; row.pop("matches", None)
            if r["conv"] and r["E"]-best <= -HIT_EV: hits.append(rep)
        return hits
    HITC = {c: h for c in INFO if INFO[c]["kept"] and all(n in D3 for n in INFO[c]["ex3"]) and (h := rank(c, INFO[c]["ex3"]))}
    extra = [(n, CH[n]["structure"], 0.01, 600, str(RUN/"d3"), DL + 600) for c in HITC for n in INFO[c]["ex"] if n not in INFO[c]["ex3"]]
    D3.update(run_stage(extra, True, "d3"))
    HITC = {c: h for c in HITC if (h := rank(c, INFO[c]["ex"]))}
    for c in INFO:                                   # mark chemistries whose D3 stage did not finish
        if INFO[c]["kept"] and not all(n in D3 for n in INFO[c]["ex3"]):
            for rep in INFO[c]["kept"]: ROWS["".join(c), rep]["status"] = "dropped_budget_D3"

    # ---- 4. decomposition AX + MX3 -> AMX4 for hit chemistries ----
    bjobs, BIN = [], {}
    for c in HITC:
        for key, form in (("AX", rform({c[0]: 1, c[2]: 1})), ("MX3", rform({c[1]: 1, c[2]: 3}))):
            es = [e for e in fetch(form) if e["nsites"] <= MAXN]
            if es:
                e = min(es, key=lambda e: e["_alexandria_hull_distance"]); BIN[c, key] = e
                bjobs.append((f"bin_{e['id']}", struct(e).as_dict(), 0.01, 600, str(RUN/"d3"), time.time() + 900))
    D3.update(run_stage(list({j[0]: j for j in bjobs}.values()), True, "d3"))
    HITS = []
    for c, reps in HITC.items():
        cn = "".join(c); ax, mx = BIN.get((c, "AX")), BIN.get((c, "MX3"))
        rb = [D3.get(f"bin_{b['id']}") if b else None for b in (ax, mx)]
        seen = []                                    # hits that relaxed into the same structure are reported once (lowest first)
        for rep in sorted(reps, key=lambda q: D3[f"tp_{cn}_{q}"]["E"]):
            row = ROWS[cn, rep]; r = D3[f"tp_{cn}_{rep}"]
            if all(ok(x) for x in rb):
                row["E_decomp_meV"] = round(1000*(r["E"]-(2*rb[0]["E"]+4*rb[1]["E"])/6), 2)
                row["binaries"] = f"{ax['id']}+{mx['id']}"; row["binaries_conv"] = bool(rb[0]["conv"] and rb[1]["conv"])
            sr = S(r); d = next((q for q, ss in seen if same(sr, ss)), None)
            if d: row.update(status="hit_duplicate_of_other_transplant", matches=d); continue
            seen.append((rep, sr)); row["status"] = "HIT_proposed_model"; HITS.append((cn, rep))
    HITS.sort(key=lambda k: (ROWS[k]["dE_D3_meV"], ROWS[k].get("E_decomp_meV", 1e9)))

    # ---- candidate folders (rebuilt each run so they always mirror the current hit list) ----
    for p in CAND.glob("TP-*"): shutil.rmtree(p)
    cnt = {}
    for cn, rep in HITS:
        row = ROWS[cn, rep]; c = next(c for c in TARGETS if "".join(c) == cn); cnt[cn] = cnt.get(cn, 0) + 1
        d = CAND/f"TP-{cn}4-{cnt[cn]:02d}"; d.mkdir(); row["candidate"] = d.name
        Structure.from_dict(D3[f"tp_{cn}_{rep}"]["structure"]).to(filename=str(d/"relaxed_mace_d3.cif"))
        _, prov = transplant(rep, c)
        json.dump({"label": "PROPOSED MODEL - prototype transplant ranked by MACE-MPA-0 + D3(BJ); not DFT, not a stability or discovery claim",
                   "formula": cn + "4", **prov, "relaxation": "mace_mp medium-mpa-0, cpu, float64, dispersion=True (D3 BJ, PBE), FIRE + FrechetCellFilter, fmax 0.01",
                   **{k: row.get(k) for k in ("dE_cheap_meV", "dE_D3_meV", "best_existing_D3", "n_existing_D3", "E_decomp_meV", "binaries", "sg_D3", "V_D3", "D3_conv")},
                   "identity": "no species-aware StructureMatcher match (default 0.2/0.3/5, loose 0.3/0.5/10) to any Alexandria entry of this composition (PBE, MACE and MACE+D3 geometries)",
                   "not_checked": "MP, OQMD, GNoME, ICSD, literature; phonons; DFT"}, open(d/"provenance.json", "w"), indent=1)

    # ---- outputs ----
    cols = ["chem", "proto", "proto_formula", "n_atoms", "n_existing", "status", "matches", "volume_scale", "dE_cheap_meV", "cheap_conv", "V_cheap",
            "dE_D3_meV", "best_existing_D3", "n_existing_D3", "D3_conv", "V_D3", "sg_D3", "E_decomp_meV", "binaries", "binaries_conv", "candidate"]
    with open(D/"results"/"transplant_cycle3.csv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, extrasaction="ignore"); w.writeheader(); w.writerows(ROWS.values())
    tf = CACHE/"walltime.json"; W = json.load(open(tf)) if tf.exists() else []
    W.append(time.time()-T0); json.dump(W, open(tf, "w"))
    st = {}
    for r in ROWS.values(): st[r["status"]] = st.get(r["status"], 0) + 1
    dropped = sorted({r["chem"] for r in ROWS.values() if r["status"].startswith("dropped")})
    e = lambda n: D3[n]["E"] if ok(D3.get(n)) else float("nan")
    ga = GOF.get("agm005743614"); ra = ROWS.get(("NaBiBr", ga), {})
    a1 = 1000*(e(f"tp_NaBiBr_{ga}") - e("ex_agm069058740")); a2 = 1000*(e("ex_agm088474086") - e("ex_agm005743614"))
    L = ["# Cycle 3: prototype-transplant scan, A M X4 (A = Li, Na; M = Al, Ga, In, Sc, Y, Sb, Bi; X = Cl, Br, I)", "",
         "Every structure here is a **PROPOSED MODEL** ranked by a machine-learned potential (MACE-MPA-0, with D3 for the final ranking). "
         "Nothing below is a stability or discovery claim; no DFT was run.", "",
         "## Counts", "",
         f"- Alexandria entries downloaded (<= {MAXN} atoms): {len(BYID)} across {len(ENT)} A'M'X'4 chemistries.",
         f"- Prototype pool (hull <= {1000*POOL_HULL:.0f} meV/atom): {len(pool)} entries -> **{len(REPS)} distinct prototypes** "
         "(role-preserving StructureMatcher at default tolerance: A', M', X' are mapped to fixed roles before matching, which is stricter than fit_anonymous).",
         f"- Chemistry x prototype pairs: {len(ROWS)}; already listed for that chemistry: {st.get('already_listed', 0)}; "
         f"**transplants attempted (cheap MACE pass): {n_att}**.",
         f"- Kept for MACE+D3 (within {1000*KEEP_EV:.0f} meV/atom of best existing in the cheap pass, not matching an existing entry, de-duplicated): {n_kept}.",
         f"- Row status counts: {', '.join(f'{k}: {v}' for k, v in sorted(st.items()))}.",
         f"- Dropped for budget (chemistries not finished before the stage deadline): {', '.join(dropped) if dropped else 'none'}.",
         f"- Prototype scope this run: {len(USE)} of {len(REPS)} prototypes (<= {MAXP} atoms{', group must contain a Cl/Br/I member' if NONF else ''}); "
         f"{st.get('skipped_prototype_scope', 0)} chemistry x prototype pairs skipped for that reason.", "",
         f"## Hits: survive identity check and are >= {1000*HIT_EV:.0f} meV/atom below the best existing entry at MACE+D3", "",
         "dE = E(transplant) - E(lowest Alexandria entry of the same composition), both MACE+D3 relaxed; for hit chemistries every existing entry "
         "(<= 24 atoms) was relaxed with D3, not only the three lowest. E_decomp = E(AMX4) - [2 E(AX) + 4 E(MX3)]/6 per atom, positive = unstable against the binaries.", "",
         "| candidate | composition | prototype source | dE no-D3 cheap (meV/atom) | dE MACE+D3 (meV/atom) | best existing (D3) | E_decomp D3 (meV/atom) | binaries AX + MX3 | SG | V (A^3/atom) |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for k in HITS:
        r = ROWS[k]
        L.append(f'| {r["candidate"]} | {r["chem"]} | {r["proto"]} ({r["proto_formula"]}) | {r["dE_cheap_meV"]:+.1f} | {r["dE_D3_meV"]:+.1f} | {r["best_existing_D3"]} | '
                 f'{r.get("E_decomp_meV", float("nan")):+.1f} | {r.get("binaries", "-")} | {r["sg_D3"]} | {r["V_D3"]} |')
    if not HITS: L.append("| none | | | | | | | | | |")
    near = sorted((r for r in ROWS.values() if r["status"] == "survives_identity"), key=lambda r: r["dE_D3_meV"])
    L += ["", f"Survivors of the identity check that are NOT hits (dE > -{1000*HIT_EV:.0f} meV/atom or unconverged): {len(near)}; lowest ten:", "",
          "| composition | prototype source | dE MACE+D3 | converged |", "|---|---|---|---|"]
    L += [f'| {r["chem"]} | {r["proto"]} ({r["proto_formula"]}) | {r["dE_D3_meV"]:+.1f} | {r["D3_conv"]} |' for r in near[:10]]
    L += ["", "## Sanity anchors", "",
          f"- (a) NaBiBr4 on the NaBiCl4/KBiBr4 prototype (pool representative {ga}) vs agm069058740, MACE+D3: **{a1:+.1f} meV/atom** "
          f"(skeptic cycle 2: +2.6). Row status: {ra.get('status')}; cheap no-D3 dE {ra.get('dE_cheap_meV')} (skeptic: -8.1 at fmax 0.01).",
          f"- (b) NaBiCl4 agm088474086 vs agm005743614, MACE+D3: **{a2:+.1f} meV/atom** (skeptic cycle 2: -12.6), i.e. agm088474086 is "
          f"{'lower' if a2 < 0 else 'NOT lower'}.", "",
          "## Wall time", "", f"- This run: {W[-1]/60:.1f} min; all runs of this script so far (resumed runs skip finished work): {sum(W)/60:.1f} min over {len(W)} run(s). {NPROC} single-thread workers.", "",
          "## Caveats", "",
          "- In-sample ML potential: MACE-MPA-0 was trained on Alexandria/MPtrj, so the existing entries are in its training set and the transplants are not. "
          "A margin of a few meV/atom between an in-sample and an out-of-sample structure is within model error.",
          "- D3 damping: Becke-Johnson with PBE parameters (mace_mp default, torch-dftd). MACE was not trained with D3; D3 shifts relative energies of these halides by 10-20 meV/atom, "
          "and zero-damping D3 (as in the existing QE inputs) was not tested.",
          "- Alexandria PBE hull distances are not comparable to MACE or MACE+D3 energies. Decomposition energies here use only AX + MX3 (single lowest-hull polymorph each); "
          "other ternaries (e.g. A3MX6) are not considered, so E_decomp is not a hull distance.",
          "- No phonons, no finite temperature, no symmetry constraints; FIRE can stall on saddle points. 'Converged' means fmax < 0.01 eV/A within 600 steps.",
          "- Triage is done without D3 (15 meV/atom window, fmax 0.05, <= 150 steps): a prototype that only becomes competitive once dispersion is on can be missed, "
          "and only the three lowest no-D3 existing entries are re-relaxed with D3 unless the chemistry has a hit.",
          "- A weak baseline makes a large dE cheap: check the Alexandria hull distance of the best existing entry (e.g. NaSbCl4 has only two listed entries, at +46 and +148 meV/atom PBE). "
          "A transplant can also relax far from its source prototype; the CIF, not the prototype name, is the proposed model.",
          "- Cells limited to <= 24 atoms; 'not listed' means not listed in Alexandria PBE for that reduced formula. MP, OQMD, GNoME, ICSD and literature were not checked.",
          "- Per-row data: `discovery/results/transplant_cycle3.csv`; raw relaxations: `discovery/results/transplant_cycle3/{cheap,d3}/`."]
    open(D/"results"/"transplant_cycle3.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L[4:12])); print("anchors", a1, a2, "hits", [(k, ROWS[k]["dE_D3_meV"], ROWS[k].get("E_decomp_meV")) for k in HITS]); print("DONE", time.time()-T0)
