"""Widen Part B/C: prototype-transplant screen (transplant_scan.py logic, generalised to A_a M_m X_x families) for
compounds that are experimentally reported but have no published structure. Everything is a PROPOSED MODEL ranked by
MACE-MPA-0 (+D3); nothing here is a stability or discovery claim. Resumable (one json per relaxation).
Run from repo root:  taskpolicy -b discovery/.venv/bin/python -I discovery/scripts/widen/widen_scan.py [cheap_min=70] [d3_min=45]
Workers: 2 while `pgrep -f nbc_search` finds a process, else 4 (re-checked at each stage)."""
import sys, json, time, subprocess
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from multiprocessing import Pool
from pymatgen.core import Structure, Element, Composition
import transplant_scan as ts                      # reuse fetch / relax workers / matchers unchanged

OLD = Path("discovery/cache"); W = Path("discovery/results/widen"); MY = W/"cache"; CAND = Path("discovery/candidates")
A_P, M_P, X_P = ts.A_P, ts.M_ALL, ts.X_P          # same prototype chemistries as the cycle-3 scan
POOL_HULL, MAXN, KEEP_EV, MAXKEEP, N_EX_D3, SHORT = 0.040, 24, 0.015, 2, 3, 12
HIT_DE, HIT_DECOMP, NEAR_HULL, MAX_COMP_N = -0.005, -0.003, 0.005, 40
TARGETS = [("K", "Bi", "Cl", (1, 1, 4)), ("Rb", "Bi", "Cl", (1, 1, 4)), ("K", "Bi", "Cl", (3, 1, 6)), ("Cs", "Bi", "Cl", (1, 1, 4)),
           ("K", "Sb", "Cl", (1, 1, 4)), ("Rb", "Bi", "Cl", (3, 1, 6)), ("K", "Bi", "Br", (1, 1, 4)), ("Rb", "Sb", "Cl", (1, 1, 4))]
NT = 3     # note: only the first NT targets and SHORT prototypes each (relaxations under taskpolicy -b ran ~8x slower than in cycle 3);
           # raise NT / SHORT (argv 3, 4) and rerun to extend, finished relaxations are reused

def shortlist(k, groups, ent, n):
    """Prototypes most likely to matter: first those already realised by a listed Bi/Sb chloride (lowest hull first, up to n//2),
    then (A M X4 only) those with the lowest cheap dE in the cycle-3 scan of Bi/Sb chlorides, then lowest source hull."""
    import csv
    real = [g for g in groups if any(ent[i]["M"] in ("Bi", "Sb") and ent[i]["X"] == "Cl" for i in g if i in ent)][:n//2]
    sc = {}
    if k == (1, 1, 4):
        for r in csv.DictReader(open("discovery/results/transplant_cycle3.csv")):
            if r["chem"] in ("NaBiCl4", "NaSbCl4", "LiBiCl4", "LiSbCl4") and r["dE_cheap_meV"]:
                sc[r["proto"]] = min(sc.get(r["proto"], 1e9), float(r["dE_cheap_meV"]))
    rest = sorted((g for g in groups if g not in real), key=lambda g: (sc.get(g[0], 1e9), ent[g[0]]["_alexandria_hull_distance"]))
    return real + rest[:n-len(real)]

def form(c, k): return ts.rform({c[0]: k[0], c[1]: k[1], c[2]: k[2]})
def name(c, k): return "".join(e + (str(n) if n > 1 else "") for e, n in zip(c, k))
def fetch(f):                                      # read-only reuse of the cycle-3 cache; new downloads go to the widen cache
    p = OLD/f"alex_{f}.json"
    if p.exists(): return json.load(open(p))
    ts.CACHE = MY; return ts.fetch(f)
def nproc(): return 2 if subprocess.run(["pgrep", "-f", "nbc_search"], capture_output=True).stdout.strip() else 4

def run_stage(jobs, d3, tag):
    out = W/tag; out.mkdir(parents=True, exist_ok=True)
    todo = [j for j in jobs if not (out/f"{j[0]}.json").exists()]; n = nproc()
    print(f"[{tag}] {len(jobs)} jobs, {len(todo)} to run, {n} workers t={time.time()-T0:.0f}s", flush=True)
    if todo:
        with Pool(n, initializer=ts._init, initargs=(d3,)) as p:
            for i, _ in enumerate(p.imap(ts._relax, todo, chunksize=1)):      # ordered: earlier targets finish first
                if i % 50 == 0: print(f"  [{tag}] {i}/{len(todo)} t={time.time()-T0:.0f}s", flush=True)
    return {j[0]: json.load(open(out/f"{j[0]}.json")) for j in jobs if (out/f"{j[0]}.json").exists()}

def line_hull(t0, phases):
    """Lowest energy per cation at cation fraction t0 from two-phase mixtures on the AX-MX3 line. phases: [(t, E_per_cation)]."""
    best = min((e for t, e in phases if abs(t-t0) < 1e-9), default=float("inf"))
    for tp, ep in phases:
        for tq, eq in phases:
            if tp < t0 - 1e-9 and tq > t0 + 1e-9: best = min(best, ep + (eq-ep)*(t0-tp)/(tq-tp))
    return best
assert abs(line_hull(0.5, [(0, 0.0), (1, 0.0), (0.25, -1.0)]) + 2/3) < 1e-9   # self-check

if __name__ == "__main__":
    T0 = time.time()
    cheap_min = float(sys.argv[1]) if len(sys.argv) > 1 else 70
    d3_min = float(sys.argv[2]) if len(sys.argv) > 2 else 45
    S = lambda r: ts.prim(Structure.from_dict(r["structure"]))
    NT = int(sys.argv[3]) if len(sys.argv) > 3 else NT; SHORT = int(sys.argv[4]) if len(sys.argv) > 4 else SHORT; TARGETS = TARGETS[:NT]
    # ---- prototype pools per stoichiometry family ----
    FAM = {}
    for k in sorted({t[3] for t in TARGETS}):
        ent = {}
        for A in A_P:
            for M in M_P:
                for X in X_P:
                    for e in fetch(form((A, M, X), k)):
                        if e["nsites"] <= MAXN: e.update(A=A, M=M, X=X); ent[e["id"]] = e
        def dummy(e):
            s = ts.struct(e); s.replace_species({e["A"]: "Na", e["M"]: "Bi", e["X"]: "Cl"}); return ts.prim(s)
        pool = sorted((e for e in ent.values() if e["_alexandria_hull_distance"] <= POOL_HULL), key=lambda e: e["_alexandria_hull_distance"])
        gf = MY/f"groups_{'-'.join(map(str, k))}.json"
        if gf.exists(): groups = json.load(open(gf))
        elif k == (1, 1, 4):   # reuse the cycle-3 grouping (pool hull <= 30 meV, full pairwise matching); regrouping at 40 meV was too slow here
            groups = json.load(open(OLD/"prototype_groups.json")); json.dump(groups, open(gf, "w"))
        else:   # note: greedy role-preserving grouping inside (nsites, Alexandria space group) buckets; same type listed under two
                # space groups stays as two prototypes (costs a duplicate relaxation, removed later by the identity checks)
            groups, dum = [], {}
            for e in pool:
                d = dum[e["id"]] = dummy(e); key = (e["nsites"], e["_alexandria_space_group"])
                g = next((g for g in groups if (ent[g[0]]["nsites"], ent[g[0]]["_alexandria_space_group"]) == key and ts.SM.fit(dum[g[0]], d)), None)
                g.append(e["id"]) if g else groups.append([e["id"]])
            json.dump(groups, open(gf, "w"))
        groups = sorted((g for g in groups if g[0] in ent), key=lambda g: ent[g[0]]["_alexandria_hull_distance"])
        n_all = len(groups); groups = shortlist(k, groups, ent, SHORT)
        FAM[k] = dict(ent=ent, groups=groups, dummy=dummy, n_all=n_all)
        print(f"family {k}: {len(ent)} entries <= {MAXN} atoms, pool {len(pool)} -> {n_all} prototypes, {len(groups)} shortlisted t={time.time()-T0:.0f}s", flush=True)

    # ---- transplants + cheap pass ----
    DL = time.time() + cheap_min*60; cheap, META = [], {}
    for c4 in TARGETS:
        c, k = c4[:3], c4[3]; tn = name(c, k); F = FAM[k]
        ex = [e for e in fetch(form(c, k)) if e["nsites"] <= 40]             # listed baselines: up to 40 atoms
        exd = {e["id"]: ts.prim(ts.struct(e)) for e in ex}
        vt = float(np.mean([ts.struct(e).volume/e["nsites"] for e in ex])) if ex else None
        tps = {}
        for g in F["groups"]:
            e = F["ent"][g[0]]
            if any(i in exd for i in g): continue                               # prototype already listed for this composition
            s = ts.struct(e); v0 = s.volume/len(s)
            if vt is None:
                r3 = lambda a, m, x: k[0]*Element(a).average_ionic_radius**3 + k[1]*Element(m).average_ionic_radius**3 + k[2]*Element(x).average_ionic_radius**3
                v = v0*float(r3(*c)/r3(e["A"], e["M"], e["X"]))
            else: v = vt
            mp = {e["A"]: c[0], e["M"]: c[1], e["X"]: c[2]}
            s.replace_species({a: b for a, b in mp.items() if a != b}); s.scale_lattice(v*len(s))
            if any(len(ts.prim(s)) == len(x) and ts.SM.fit(ts.prim(s), x) for x in exd.values()): continue
            tps[g[0]] = dict(source_id=g[0], source_formula=name((e["A"], e["M"], e["X"]), k), source_hull_eV=e["_alexandria_hull_distance"],
                             source_sg=e["_alexandria_space_group"], substitution=mp, volume_scale=v/v0, prototype_group=g)
            cheap.append((f"tp_{tn}_{g[0]}", s.as_dict(), 0.05, 150, str(W/"cheap"), DL))
        for e in ex: cheap.append((f"ex_{e['id']}", ts.struct(e).as_dict(), 0.05, 150, str(W/"cheap"), DL))
        META[tn] = dict(c=c, k=k, ex=ex, exd=exd, tps=tps)
        print(f"{tn}: {len(ex)} listed, {len(tps)} transplants", flush=True)
    CH = run_stage(list({j[0]: j for j in cheap}.values()), False, "cheap")

    # ---- triage -> D3 (survivors + lowest listed after a 0.03 A rattle) ----
    DL = time.time() + d3_min*60; d3 = []
    for tn, m in META.items():
        exn = [f"ex_{e['id']}" for e in m["ex"]]; tpn = [f"tp_{tn}_{r}" for r in m["tps"]]
        m["complete_cheap"] = all(n in CH for n in exn + tpn)
        exok = sorted((n for n in exn if ts.ok(CH.get(n))), key=lambda n: CH[n]["E"])
        if not exok: m["status"] = "no listed reference relaxed"; continue
        best = CH[exok[0]]["E"]; exS = [x for n in exok for x in (m["exd"][n[3:]], S(CH[n]))]; kept = []
        for n in sorted((n for n in tpn if ts.ok(CH.get(n))), key=lambda n: CH[n]["E"]):
            r = CH[n]
            if r["E"]-best > KEEP_EV: break
            s = S(r)
            if any(ts.same(s, x) for x in exS) or any(ts.same(s, x) == "default" for _, x in kept): continue
            kept.append((n, s))
        m.update(exok=exok, kept=[n for n, _ in kept][:MAXKEEP], n_kept_all=len(kept), ex3=exok[:N_EX_D3])
        for n in m["ex3"]:
            np.random.seed(0); s = Structure.from_dict(CH[n]["structure"]); s.perturb(0.03)
            d3.append((n + "_rattled", s.as_dict(), 0.01, 600, str(W/"d3"), DL))
        for n in m["kept"]: d3.append((n, CH[n]["structure"], 0.01, 600, str(W/"d3"), DL))
    D3 = run_stage(list({j[0]: j for j in d3}.values()), True, "d3")

    # ---- competitors on the AX - MX3 line: binaries (2 lowest-hull polymorphs each) + listed ternaries within 5 meV of the hull ----
    comp, COMP = [], {}
    for tn, m in META.items():
        if not m.get("kept"): continue
        c, k = m["c"], m["k"]; ph = []
        for f, t in ((ts.rform({c[0]: 1, c[2]: 1}), 0.0), (ts.rform({c[1]: 1, c[2]: 3}), 1.0)):
            for e in sorted((e for e in fetch(f) if e["nsites"] <= MAXN), key=lambda e: e["_alexandria_hull_distance"])[:1]: ph.append((f, t, e))
        seen = {}
        for e in json.load(open(MY/f"alex_{c[0]}-{c[1]}-{c[2]}.json")):
            d = Composition(e["chemical_formula_reduced"]).get_el_amt_dict(); nA, nM, nX = d.get(c[0], 0), d.get(c[1], 0), d.get(c[2], 0)
            if abs(nX - nA - 3*nM) > 1e-6 or e["_alexandria_hull_distance"] > NEAR_HULL or abs(nM/(nA+nM) - k[1]/(k[0]+k[1])) < 1e-6: continue
            f = e["chemical_formula_reduced"]
            if f not in seen or e["_alexandria_hull_distance"] < seen[f][1]["_alexandria_hull_distance"]: seen[f] = (nM/(nA+nM), e)
        m["skipped_competitors"] = []
        for f, (t, e) in seen.items():
            if e["nsites"] > MAX_COMP_N: m["skipped_competitors"].append(f"{f} ({e['id']}, {e['nsites']} atoms)"); continue
            full = next(x for x in fetch(f) if x["id"] == e["id"]); ph.append((f, t, full))
        COMP[tn] = ph
        for f, t, e in ph:
            if not ts.ok(D3.get(f"ex_{e['id']}_rattled")): comp.append((f"cp_{e['id']}", ts.struct(e).as_dict(), 0.01, 600, str(W/"d3"), time.time() + 1200))
    D3.update(run_stage(list({j[0]: j for j in comp}.values()), True, "d3"))
    for tn in COMP:
        for f, t, e in COMP[tn]:
            if ts.ok(D3.get(f"ex_{e['id']}_rattled")): D3[f"cp_{e['id']}"] = D3[f"ex_{e['id']}_rattled"]

    # ---- ranking, identity, decomposition, hand-off ----
    OUT, L, cnt = {}, [], {}
    BEST = {}   # lowest MACE+D3 energy found for each screened composition (listed or transplant), so a competitor is never weaker than our own best model of it
    for tn, m in META.items():
        es = [D3[n]["E"] for n in [x + "_rattled" for x in m.get("ex3", [])] + m.get("kept", []) if ts.ok(D3.get(n)) and D3[n]["conv"]]
        if es: BEST[Composition(tn).reduced_formula] = min(es)
    for tn, m in META.items():
        c, k = m["c"], m["k"]; nat = sum(k); ncat = k[0]+k[1]; o = OUT[tn] = dict(n_listed=len(m["ex"]), n_transplants=len(m["tps"]),
            complete_cheap=m.get("complete_cheap"), n_kept=m.get("n_kept_all"), status=m.get("status", ""), skipped_competitors=m.get("skipped_competitors", []))
        exr = [n for n in m.get("ex3", []) if ts.ok(D3.get(n + "_rattled"))]
        if not exr: o["status"] = o["status"] or "D3 of listed entries not finished"; continue
        bn = min(exr, key=lambda n: D3[n + "_rattled"]["E"]); rb = D3[bn + "_rattled"]; sb = S(rb)
        o["best_listed"] = dict(id=bn[3:], E=rb["E"], V=rb["V"], sg=ts.sgs(sb), n=rb["n"], conv=rb["conv"],
                                alex_hull_meV=1000*next(e["_alexandria_hull_distance"] for e in m["ex"] if e["id"] == bn[3:]),
                                listed_D3_meV={n[3:]: round(1000*(D3[n + "_rattled"]["E"]-rb["E"]), 2) for n in exr})
        ph = [(t, min(D3[f"cp_{e['id']}"]["E"], BEST.get(Composition(f).reduced_formula, 9e9))*sum(Composition(f).get_el_amt_dict().values())/(Composition(f)[c[0]]+Composition(f)[c[1]]), f, e["id"])
              for f, t, e in COMP.get(tn, []) if ts.ok(D3.get(f"cp_{e['id']}"))]
        o["competitors"] = [dict(formula=f, id=i, conv=D3[f"cp_{i}"]["conv"], n=D3[f"cp_{i}"]["n"]) for _, _, f, i in ph]
        t0 = k[1]/ncat
        def decomp(E, binaries_only):
            p = [(t, e) for t, e, f, i in ph if not binaries_only or t in (0.0, 1.0)]
            if not {0.0, 1.0} <= {t for t, _ in p}: return None
            return 1000*(E*nat/ncat - line_hull(t0, [x for x in p if abs(x[0]-t0) > 1e-9]))*ncat/nat
        o["best_listed"].update(decomp_binaries_meV=decomp(rb["E"], True), decomp_all_meV=decomp(rb["E"], False))
        exS = [x for n in m["exok"] for x in [m["exd"][n[3:]], S(CH[n])] + ([S(D3[n + "_rattled"])] if ts.ok(D3.get(n + "_rattled")) else [])]
        rows = []
        for n in m["kept"]:
            r = D3.get(n)
            if not ts.ok(r): rows.append(dict(name=n, status="D3 not finished/failed")); continue
            s = S(r); idm = next((t for x in exS if (t := ts.same(s, x))), "")
            rows.append(dict(name=n, proto=m["tps"][n.split("_")[-1]]["source_formula"], dE_meV=1000*(r["E"]-rb["E"]), dE_cheap_meV=1000*(CH[n]["E"]-CH[m["exok"][0]]["E"]),
                             conv=r["conv"], V=r["V"], sg=ts.sgs(s), n=r["n"], matches_listed=idm, decomp_binaries_meV=decomp(r["E"], True), decomp_all_meV=decomp(r["E"], False)))
        rows.sort(key=lambda r: r.get("dE_meV", 9e9)); o["transplants"] = rows; seen = []
        for r in rows:
            if "dE_meV" not in r or r["matches_listed"] or not r["conv"]: continue
            s = S(D3[r["name"]])
            if any(ts.same(s, x) for x in seen): r["duplicate"] = True; continue
            seen.append(s)
            if r["dE_meV"] <= 1000*HIT_DE and r["decomp_all_meV"] is not None and r["decomp_all_meV"] <= 1000*HIT_DECOMP:
                cnt[tn] = cnt.get(tn, 0) + 1; d = CAND/f"W-{tn}-{cnt[tn]:02d}"; d.mkdir(parents=True, exist_ok=True); r["candidate"] = d.name
                Structure.from_dict(D3[r["name"]]["structure"]).to(filename=str(d/"relaxed_mace_d3.cif"))
                json.dump({"label": "PROPOSED MODEL - prototype transplant ranked by MACE-MPA-0 + D3(BJ); not DFT, not a stability or discovery claim",
                           "formula": tn, **m["tps"][r["name"].split("_")[-1]],
                           "relaxation": "mace_mp medium-mpa-0, cpu, float64, dispersion=True (D3 BJ, PBE), FIRE + FrechetCellFilter, fmax 0.01, <= 600 steps",
                           **{q: r[q] for q in ("dE_meV", "dE_cheap_meV", "decomp_binaries_meV", "decomp_all_meV", "sg", "V", "n", "conv")},
                           "best_listed": o["best_listed"], "competitors_relaxed": o["competitors"], "competitors_not_relaxed": o["skipped_competitors"],
                           "identity": "no species-aware StructureMatcher match (default 0.2/0.3/5, loose 0.3/0.5/10) to any listed Alexandria entry (<= 40 atoms) of this composition (PBE, MACE, MACE+D3 geometries)",
                           "not_checked": "ICSD, Pauling File, COD, OQMD, GNoME; phonons; DFT"}, open(d/"provenance.json", "w"), indent=1)
    json.dump(OUT, open(W/"scan_summary.json", "w"), indent=1, default=str)
    f1 = lambda v: "n/a" if v is None else f"{v:+.1f}"
    for tn, o in OUT.items():
        L.append(f"### {tn}: {o['n_listed']} listed, {o['n_transplants']} transplants, cheap pass complete: {o['complete_cheap']}, kept {o['n_kept']} {o['status']}")
        if "best_listed" not in o: L.append(""); continue
        b = o["best_listed"]
        L += [f"- best listed (MACE+D3, rattled): {b['id']} sg {b['sg']} n={b['n']} V={b['V']:.2f} A^3/atom, Alexandria hull {b['alex_hull_meV']:.1f} meV; "
              f"decomposition vs binaries {f1(b['decomp_binaries_meV'])}, vs all competitors {f1(b['decomp_all_meV'])} meV/atom; other listed: {b['listed_D3_meV']}",
              f"- competitors relaxed: {', '.join(x['formula']+':'+x['id']+('' if x['conv'] else '(unconv)') for x in o['competitors'])}; not relaxed (too large): {o['skipped_competitors'] or 'none'}", "",
              "| transplant | prototype | dE vs best listed | dE cheap | vs binaries | vs all competitors | conv | matches listed | SG | n | V | candidate |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in o["transplants"]:
            L.append(f"| {r['name']} | {r.get('proto','')} | {f1(r.get('dE_meV'))} | {f1(r.get('dE_cheap_meV'))} | {f1(r.get('decomp_binaries_meV'))} | {f1(r.get('decomp_all_meV'))} | {r.get('conv')} | "
                     f"{r.get('matches_listed') or ('dup' if r.get('duplicate') else 'no')} | {r.get('sg','')} | {r.get('n','')} | {r.get('V',0):.2f} | {r.get('candidate','')} |")
        L.append("")
    open(W/"partB.md", "w").write("\n".join(L) + "\n"); print("\n".join(L)); print("DONE", time.time()-T0)
