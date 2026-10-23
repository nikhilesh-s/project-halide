"""Writes discovery/results/nbc_search.md, discovery/candidates/NBC-nn/ (structures at or below the winner) and QE inputs for
the top <= 4 new contenders in overnight/later/nbc_rivals/. Reads caches only (plus the reference relaxations if missing)."""
import sys, os, glob, collections, shutil; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import s2_airss as A2, s3_perturb as P3
from d3_final import WIN, load, contenders
import make_inputs as mi
def cn(s, el, cut):
    out = collections.Counter()
    for site in s:
        if site.specie.symbol == el:
            d = sorted(x.nn_distance for x in s.get_neighbors(site, cut) if x.specie.symbol == "Cl")
            out[f"CN{len(d)} ({d[0]:.2f}-{d[-1]:.2f})" if d else "CN0"] += 1
    return "; ".join(out)
base = lambda m: m.split("(")[0].replace("_sk3geom", "").replace("_MLgeom", "").replace("_DFTgeom", "")

if __name__ == "__main__":
    ch, d3 = ref_results(); w = d3[WIN]["E"]; wc = ch[WIN]["E"]; lib = ref_lib(ch, d3)
    PR = json.load(open(C/"protos.json")); PID = {p["id"]: p for p in PR["protos"]}; G = json.load(open(C/"airss_starts.json"))
    sel = json.load(open(C/"selection.json")) if (C/"selection.json").exists() else dict(k=0, k_rms=0, rows=[], n_pool=0, n_window=0)
    s1, s2, s3, cd = load("s1_cheap"), load("airss_cheap"), load("s3_cheap"), {n: r for n, r in load("cand_d3").items() if ok(r)}
    def prov(n):
        n = n[2:] if n.startswith("c_") else n
        if n.startswith("tp_"): p = PID[n[3:]]; return f"S1 transplant of {p['form']} {p['id']} ({p['src']}, {p['n']} atoms, hull {1000*p['hull']:.0f} meV)"
        if n.startswith("rs_"): return f"S2 random start {n} (Z={G[n]['Z']}, start symmetry {G[n]['sg_start']})"
        return f"S3 perturbation {n}"
    def redisc(res):   # which known structures were re-found by the cheap relaxations
        c = collections.Counter()
        for n, r in res.items():
            if ok(r) and r["E"]-wc < 0.06 and (m := identify(S(r), lib)): c[base(m)] += 1
        return c
    ok1 = {n: r for n, r in s1.items() if ok(r)}; ok2 = {n: r for n, r in s2.items() if ok(r)}
    r1, r2 = redisc(ok1), redisc(ok2); cl2 = A2.distinct(ok2); KNOWN = sorted({base(k) for k in lib})
    def ph(n):
        f = C/"ph"/f"{n}.json"
        if not f.exists(): return "not run / not finished"
        p = json.load(open(f)); return f"min {p['min_THz']:+.2f} THz at q={p['q_min']} ({p['n_below_m0p3']} modes < -0.3)"
    L = ["# NaBiCl4: adversarial search for anything below the KBiF4-type winner (cycle 5)", "",
         "Everything here is a **PROPOSED MODEL** ranked by a machine-learned potential (MACE-MPA-0 medium, + D3(BJ) float64 for the final numbers). No DFT was run in this cycle; nothing is a stability or discovery claim.", "",
         "Reference = the DFT-tested winner (KBiF4-type transplant of agm077914577, DFT-relaxed geometry from `overnight/night2/jobs/01_NaBiCl4_TP01_KBiF4type_vcrelax.out`, re-relaxed with MACE+D3 after a 0.03 A rattle). All dE are meV/atom relative to it.", "",
         "## Calibration (how far to trust the ML ranking)", "",
         f"- MACE+D3: winner - agm088474086 = **{1000*(w-d3['agm088474086']['E']):+.2f} meV/atom**; DFT (E31, 80 Ry): **-5.71**.",
         f"- MACE+D3: winner - agm084661913 (its Fdd2 minimum) = {1000*(w-min(d3[k]['E'] for k in d3 if k.startswith('agm084661913'))):+.2f}; DFT: -9.56. MACE+D3: winner - agm005743614 (C2/c) = {1000*(w-d3['agm005743614']['E']):+.2f}; DFT: -25.84.",
         "- So the model reproduces one gap to 0.2 meV/atom and misses the two others by about 7.5; skeptic cycle 3 measured an RMS of 3.4 meV/atom for Na-Bi-Cl. **Any ML margin below roughly 5-8 meV/atom is not a ranking**; it is a reason to run DFT.", "",
         "## Reference set at the same level", "", "| structure | atoms | dE MACE+D3 | V A^3/atom | SG (0.1) | dE no-D3 cheap | Na-Cl (<= 3.3 A) | Bi-Cl (<= 3.2 A) |", "|---|---|---|---|---|---|---|---|"]
    for k in sorted(d3, key=lambda k: d3[k]["E"]):
        s = S(d3[k]); L.append(f"| {k} | {d3[k]['n']} | {1000*(d3[k]['E']-w):+.2f} | {d3[k]['V']:.2f} | {sg(s)} | {1000*(ch[k]['E']-wc):+.2f} | {cn(s, 'Na', 3.3)} | {cn(s, 'Bi', 3.2)} |")
    nb = lambda d: sum(1 for r in d.values() if "error" in r or not ok(r))
    srcs = collections.Counter(PID[n[3:]]["src"] for n in s1); nat = collections.Counter(PID[n[3:]]["n"] for n in s1)
    L += ["", "The KGaI4-type rival (24 atoms, already queued for DFT) is **below** the winner at this level; agm084661913 has two minima depending on the rattle (the lower one is the skeptic-cycle-3 geometry).", "",
          "## Throughput caveat", "",
          "`taskpolicy -b` (required: the laptop was running DFT and in use) puts the workers on throttled efficiency cores: a relaxation step cost about 10x what it did in the cycle-3 scan. The counts below are therefore far smaller than the few hundred to 1500 random starts that were asked for. The scripts are resumable: rerunning them with longer deadlines extends every stage.", "",
          "## Strategy 1: wider prototype pool", "",
          f"- Sources: Alexandria PBE (one anonymous-formula A4BC query + the cycle-3 cache) and Materials Project; ABX4 with A in Li/Na/K/Rb/Cs/Ag/Tl/Cu, B in Al/Ga/In/Tl/Sc/Y/La-Lu/Bi/Sb/As/Fe/Cr/Au, X in F/Cl/Br/I, <= 48 atoms, hull <= 60 meV/atom: {PR['counts']} entries. OQMD: {PR['oqmd']}.",
          f"- After role-preserving dedup (StructureMatcher default, inside equal primitive size + space group buckets): {PR['dup_of_old']} duplicates of the 228 prototypes already tried, {PR['dup_of_new']} duplicates among the rest, **{len(PID)} new prototypes** (by primitive atoms: {dict(sorted(collections.Counter(p['n'] for p in PID.values()).items()))}).",
          f"- Cheap pass (MACE float32, no D3, LBFGS, fmax 0.05, <= 150 steps, start volume 28 A^3/atom), <= 24 atoms and low hull first until the deadline: **{len(s1)} of {len(PID)} relaxed** ({dict(srcs)}; by atoms {dict(sorted(nat.items()))}); failed/collapsed {nb(s1)}. **{len(PID)-len(s1)} new prototypes were NOT tried** (all larger cells / higher hull).",
          f"- Re-found known NaBiCl4 structures: {dict(r1) or 'none'}.", "",
          "## Strategy 2: random structure search (AIRSS-style)", "",
          f"- Starts generated: {len(G)} (Z: {dict(sorted(collections.Counter(g['Z'] for g in G.values()).items()))}; {sum(g['sg_start'] != 'P1' for g in G.values())} with imposed symmetry; {sum(g['constraint_scale'] < 1 for g in G.values())} needed the distance constraints relaxed by 5 %). Z = 3 never came out of the generator; Z = 4 is under-sampled.",
          f"- Relaxed (MACE float32, no D3, LBFGS, fmax 0.05, <= 250 steps): **{len(s2)}** (Z: {dict(sorted(collections.Counter(G[n]['Z'] for n in s2).items()))}); usable {len(ok2)}, converged {sum(bool(r['conv']) for r in ok2.values())}.",
          f"- Distinct minima among the usable ones (same primitive size, |dE| < 4 meV/atom, StructureMatcher): **{len(cl2)}**; {sum(len(c['members']) > 1 for c in cl2)} of them were hit more than once.",
          f"- **Rediscovery rate**: winner re-found by {r2.get('WIN_KBiF4type', 0)} of {len(ok2)} starts; known structures re-found at least once: {sum(1 for k in KNOWN if r2.get(k))} of {len(KNOWN)} ({dict(r2) or 'none'}). Known set = winner, KGaI4-type and NaGaF4-type rivals, five Alexandria entries.",
          "- Reading: a search that does not re-find the known low-energy structures, and whose minima are almost all hit once, has not sampled this landscape; its failure to find something lower is weak evidence." if r2.get("WIN_KBiF4type", 0) == 0 else "- Reading: the winner was re-found, which is some evidence of sampling at this cell size; see the multiplicity above for how thin it is.", "",
          "## Strategy 3: local perturbations (supercell + rattle + strain)", "", "| job | parent | atoms | level | dE vs parent (meV/atom) | converged | outcome |", "|---|---|---|---|---|---|---|"]
    for n, k, sc, sd, st, direct in P3.plan():
        r = (cd if direct else s3).get(n); pe = (d3 if direct else ch)[k]["E"]
        if r is None or not ok(r): L.append(f"| {n} | {k} | - | {'MACE+D3' if direct else 'cheap'} | - | - | not run (deadline) or failed |"); continue
        m = identify(S(r), {k: lib[k]})
        L.append(f"| {n} | {k} | {r['n']} | {'MACE+D3' if direct else 'cheap no-D3'} | {1000*(r['E']-pe):+.2f} | {r['conv']} | {'returns to parent' if m else 'DISTINCT from parent (' + (identify(S(r), lib) or 'no known match') + ')'} |")
    low = {n: r for n, r in cd.items() if r["E"]-w <= 0.015}; rows = A2.distinct(low, tol=0.002)
    L += ["", "## MACE+D3 stage", "",
          f"- Triage score = dE_cheap + k (V_cheap - V_cheap,winner), k = {1000*sel['k']:.1f} meV/A^3 fitted on the reference set (rms {1000*sel['k_rms']:.1f} meV/atom). Reason: D3 favours dense cells; the no-D3 energy alone ranks the KGaI4-type rival at +10.7.",
          f"- Cheap pool {sel['n_pool']}; inside the 25 meV/atom score window (widened from 15 because so few candidates existed) {sel['n_window']}; distinct {len(sel['rows'])}; already-known structures {sum(bool(r['match']) for r in sel['rows'])}; sent to MACE+D3 (float64, LBFGS + cell, fmax 0.01, 0.03 A rattle): {sum(not r['match'] for r in sel['rows'])}; **finished before the deadline: {sum(1 for n in cd if n.startswith('c_'))}**.", "",
          "### Distinct structures within +15 meV/atom of the winner (MACE+D3), lowest first", "",
          "| # | dE meV/atom | atoms | V | SG 0.01 / 0.1 | identity vs winner, rivals, 5 Alexandria | Na-Cl | Bi-Cl | conv | phonon 2x2x2 | provenance |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, c in enumerate(rows):
        r = cd[c["rep"]]; s = S(r)
        L.append(f"| {i+1} | {1000*(r['E']-w):+.2f} | {len(c['p'])} | {r['V']:.2f} | {sg(s, 0.01)} / {sg(s)} | {identify(s, lib) or 'no match (new)'} | {cn(s, 'Na', 3.3)} | {cn(s, 'Bi', 3.2)} | {r['conv']} | {ph(c['rep'])} | {prov(c['rep'])}{' (+%d more)' % (len(c['members'])-1) if len(c['members']) > 1 else ''} |")
    if not rows: L.append("| none | | | | | | | | | | |")
    # ---- deliverables: candidate folders (at or below the winner) and QE inputs (top <= 4 new within +5 meV/atom) ----
    con = [c for c in contenders(0.005) if not c["match"]]; below = [c for c in con if c["dE"] <= 0]
    for p in glob.glob("discovery/candidates/NBC-[0-9]*"): shutil.rmtree(p)
    mi.OUT = mi.ROOT/"overnight/later/nbc_rivals"; mi.OUT.mkdir(parents=True, exist_ok=True)
    for f in mi.OUT.glob("*.in"): f.unlink()
    L += ["", "## Deliverables", ""]
    for i, c in enumerate(con):
        tag = f"NBC-{i+1:02d}"; s = S(cd[c["name"]]); line = f"- {tag}: dE {1000*c['dE']:+.2f} meV/atom, {c['n']} atoms, V {c['V']:.2f}, {prov(c['name'])}"
        if c["dE"] <= 0:
            d = Path(f"discovery/candidates/{tag}"); d.mkdir(); s.to(filename=str(d/"relaxed_mace_d3.cif"))
            json.dump({"label": "PROPOSED MODEL - ranked by MACE-MPA-0 + D3(BJ); not DFT, not a stability or discovery claim", "formula": "NaBiCl4", "provenance": prov(c["name"]),
                       "job": c["name"], "dE_vs_KBiF4type_winner_meV_per_atom_MACE_D3": round(1000*c["dE"], 2), "atoms_primitive": c["n"], "V_per_atom": round(c["V"], 2), "sg_0.1": sg(s),
                       "converged": c["conv"], "relaxation": "mace_mp medium-mpa-0, cpu, float64, dispersion=True (D3 BJ, PBE), LBFGS + FrechetCellFilter, fmax 0.01, from a 0.03 A rattle",
                       "identity": "no StructureMatcher match (default or loose) to the winner, the KGaI4-/NaGaF4-type rivals or the five Alexandria NaBiCl4 entries",
                       "phonon": ph(c["name"]), "model_error_note": "MACE vs PBE-D3 gaps for NaBiCl4 polymorphs differ by up to ~7.5 meV/atom; this margin is not a ranking"}, open(d/"provenance.json", "w"), indent=1)
            line += f"; folder `discovery/candidates/{tag}/`"
        if i < 4:
            mi.write(f"nbc{i+1:02d}_NaBiCl4_{c['name'][2:] if c['name'].startswith('c_') else c['name']}_vcrelax", s, "vc-relax", ecut=60.0); line += f"; QE input in `overnight/later/nbc_rivals/`"
        L.append(line)
    if not con: L.append("- No new distinct structure within +5 meV/atom of the winner at MACE+D3: no NBC folders, no QE inputs written.")
    kg = 1000*(d3["RIV_KGaI4type"]["E"]-w)
    L += ["", f"New distinct structures at or below the winner: **{len(below)}**; within +5 meV/atom: {len(con)}. Already known and still untested by DFT: KGaI4-type rival at {kg:+.2f} (24 atoms; `discovery/candidates/TP-NaBiCl4-01/`), NaGaF4-type at {1000*(d3['RIV_NaGaF4type']['E']-w):+.2f}.", "",
          "## Judgement", "",
          "- Nothing new was found below the winner, but the search was thin (see rediscovery rate and the untried prototypes), so this is weak support, not confirmation.",
          "- The real threats are the dense, Na-CN7 structures (V about 26 A^3/atom): the KGaI4-type rival (below the winner at MACE+D3) and any NBC entry listed above. Their ranking leans on D3(BJ); the QE runs use zero-damping D3, which was not tested here and may shift dense cells by several meV/atom either way.",
          "- With margins of 0-3 meV/atom against a model error of 3-8 meV/atom, 'lowest-energy NaBiCl4 structure found' is undecided until those DFT jobs run. Defensible wording today: 'lowest of the structures tested by DFT'.", "",
          "## What was NOT covered", "",
          "- Z > 4 (more than 24 atoms) in the random search; Z = 3 not sampled at all; Z = 4 thinly. Random starts: tens, not hundreds.",
          f"- {len(PID)-len(s1)} of the {len(PID)} new prototypes (mostly 18-48 atoms) never got a cheap relaxation; OQMD was unreachable (HTTP 502); ICSD, GNoME, AFLOW not queried. Prototype chemistries outside the A(+1) B(+3) X4 lists (e.g. A2+ B2+ X4, oxides/chalcogenides ABX4) were not used as templates.",
          "- Non-stoichiometric, disordered, partially occupied or defect models; other Na-Bi-Cl compositions (the hull question is separate: Na3BiCl6 + BiCl3 is covered by DFT in E28/E31).",
          "- Bi lone-pair physics beyond what MACE-MPA-0 learned from PBE; no spin-orbit coupling; no finite temperature / vibrational free energy; no pressure.",
          "- Triage was done without D3 (with a volume correction): a structure whose basin only exists with dispersion can be missed. Cheap relaxations were capped (150 / 250 steps), unconverged ones were kept in the ranking only by their last energy.",
          "- Phonons: only the structures listed with a phonon entry; 24-atom contenders (192-atom supercells) were not attempted.",
          "- Everything is in-sample-biased ML: the Alexandria entries are in MACE's training set, the new models are not.", "",
          "## Compute", "", f"- Core time in relaxations: cheap {sum(r.get('wall', 0) for d in (s1, s2, s3) for r in d.values())/60:.0f} core-min, MACE+D3 {sum(r.get('wall', 0) for r in load('cand_d3').values())/60:.0f} core-min; 4 single-thread workers under `taskpolicy -b`.",
          "- Scripts: `discovery/scripts/nbc_search/` (common, s1_fetch, s1_run, s2_airss, s3_perturb, d3_final, phonon, report, chain.sh). Caches: `discovery/cache/nbc_search/`."]
    open("discovery/results/nbc_search.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
