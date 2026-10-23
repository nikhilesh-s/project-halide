"""MACE+D3 float64 (fmax 0.01) for the distinct, not-already-known low minima of strategies 1-3.
Triage score = dE_cheap + k*(V_cheap - V_cheap(winner)), k fitted on the reference set: D3 favours dense cells, and the
no-D3 pass alone would have discarded the KGaI4-type rival (+10.7 cheap, -2.6 with D3).  argv: minutes | 'dry'"""
import sys, os, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import s2_airss as A2, s3_perturb as P3
WIN = "WIN_KBiF4type_DFTgeom"
def load(tag):
    out = {}
    for f in glob.glob(str(C/tag/"*.json")):
        r = json.load(open(f))
        if not r["name"].startswith("ref_"): out[r["name"]] = r
    return out
def pool(): return {n: r for t in ("s1_cheap", "airss_cheap", "s3_cheap") for n, r in load(t).items() if ok(r)}
def kfit(ch, d3):
    ks = [k for k in ch if k in d3 and k != WIN and ok(ch[k]) and ok(d3[k])]
    x = np.array([ch[k]["V"]-ch[WIN]["V"] for k in ks]); y = np.array([(d3[k]["E"]-d3[WIN]["E"])-(ch[k]["E"]-ch[WIN]["E"]) for k in ks])
    return float(x@y/(x@x)), float(np.sqrt(np.mean((y-x*(x@y/(x@x)))**2)))
def select(ch, d3, window=0.025):
    P = pool(); k, rms = kfit(ch, d3); w = ch[WIN]
    sc = {n: (r["E"]-w["E"]) + k*(r["V"]-w["V"]) for n, r in P.items()}
    cand = {n: P[n] for n in P if sc[n] <= window or P[n]["E"]-w["E"] <= 0.005}
    lib = ref_lib(ch, d3); rows = []
    for c in A2.distinct(cand):
        r = P[c["rep"]]; rows.append(dict(rep=c["rep"], members=c["members"], match=identify(S(r), lib), score=sc[c["rep"]], dE_cheap=r["E"]-w["E"], V=r["V"], n=len(c["p"])))
    return dict(k=k, k_rms=rms, n_pool=len(P), n_window=len(cand), rows=sorted(rows, key=lambda r: r["score"]))
def contenders(lim=0.005):
    """Distinct MACE+D3 minima within lim of the winner: list of dicts (lowest first), refs excluded from the D3 pool by construction."""
    ch, d3 = ref_results(); w = d3[WIN]["E"]; lib = ref_lib(ch, d3); R = {n: r for n, r in load("cand_d3").items() if ok(r)}
    out = []
    for c in A2.distinct({n: r for n, r in R.items() if r["E"]-w <= lim}, tol=0.002):
        r = R[c["rep"]]; out.append(dict(name=c["rep"], members=c["members"], dE=r["E"]-w, n=len(c["p"]), V=r["V"], conv=r["conv"], match=identify(S(r), lib), p=c["p"]))
    return out
if __name__ == "__main__":
    ch, d3 = ref_results(); sel = select(ch, d3); json.dump(sel, open(C/"selection.json", "w"), indent=1); P = pool()
    new = [r for r in sel["rows"] if not r["match"]]
    print(f"k={1000*sel['k']:.1f} meV/A^3 (rms {1000*sel['k_rms']:.1f}); pool {sel['n_pool']}, in window {sel['n_window']}, distinct {len(sel['rows'])}, already known {len(sel['rows'])-len(new)}, to D3 {len(new)}", flush=True)
    for r in sel["rows"][:40]: print(f"  {r['rep']:28s} n={r['n']:2d} score {1000*r['score']:+6.1f} cheap {1000*r['dE_cheap']:+6.1f} V {r['V']:.2f} x{len(r['members'])} {r['match']}")
    if sys.argv[1] == "dry": sys.exit()
    dl = time.time() + float(sys.argv[1])*60
    jobs = [J("c_" + r["rep"], rattle(prim(S(P[r["rep"]])), 0.03, 3), True, "cand_d3", dl, steps=300) for r in new]
    s3 = [(j[0], j[1], 0.01, 300, str(C/"cand_d3"), dl) for j in P3.jobs(True, dl)]   # direct D3 on rattled winner supercells
    R = run(jobs[:12] + s3[:1] + jobs[12:] + s3[1:], True, "cand_d3", _relax2)
    print("D3 finished", len(R), "of", len(jobs) + len(s3))
