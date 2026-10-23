"""Strategy 1, step 1: wider ABX4 prototype pool (Alexandria anonymous-formula query, Materials Project, OQMD if reachable,
plus the >24-atom / 30-60 meV entries already cached by transplant_scan.py). Role-preserving dedup against the 228 prototypes
already tried. Writes discovery/cache/nbc_search/protos.json (NaBiCl4 transplants at 28 A^3/atom). No API key is stored."""
import json, os, sys, time, glob, urllib.request, urllib.parse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # python -I does not add the script dir
from common import *
from pymatgen.core import Lattice
AL = ["Li", "Na", "K", "Rb", "Cs", "Ag", "Tl", "Cu"]
BL = ["Al", "Ga", "In", "Tl", "Sc", "Y", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Bi", "Sb", "As", "Fe", "Cr", "Au"]
XL = ["F", "Cl", "Br", "I"]; HULL, MAXN = 0.060, 48

def get(url, hdr={}):
    for k in range(4):
        try: return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=90))
        except Exception as e: err = e; time.sleep(3*(k+1))
    raise err

def alex():
    f = C/"alex_A4BC.json"
    if f.exists(): return json.load(open(f))
    q = lambda l: ",".join(f'"{x}"' for x in l)
    url = "https://alexandria.icams.rub.de/pbe/v1/structures?" + urllib.parse.urlencode({"page_limit": 100,
        "filter": f'chemical_formula_anonymous="A4BC" AND elements HAS ANY {q(XL)} AND elements HAS ANY {q(AL)} AND nsites<={MAXN} AND _alexandria_hull_distance<={HULL}',
        "response_fields": "lattice_vectors,cartesian_site_positions,species_at_sites,nsites,_alexandria_hull_distance"}); out = []
    while url:
        r = get(url); out += [{"id": d["id"], **d["attributes"]} for d in r["data"]]
        nx = (r.get("links") or {}).get("next"); url = nx.get("href") if isinstance(nx, dict) else nx; time.sleep(0.5)
        print("alex", len(out), flush=True)
    json.dump(out, open(f, "w")); return out

def mp():
    f = C/"mp_ABX4.json"
    if f.exists(): return json.load(open(f))
    key = open(os.path.expanduser("~/.mp_api_key")).read().strip(); out = []
    for X in XL:
        r = get("https://api.materialsproject.org/materials/summary/?" + urllib.parse.urlencode({"formula": f"*1*1{X}4", "nsites_max": MAXN,
            "energy_above_hull_max": HULL, "_fields": "material_id,formula_pretty,nsites,energy_above_hull,structure,theoretical", "_limit": 1000}),
            {"X-API-KEY": key, "accept": "application/json", "User-Agent": "curl/8"})
        out += r["data"]; print("mp", X, len(r["data"]), flush=True)
    json.dump(out, open(f, "w")); return out

def roles(s):
    c = s.composition.reduced_composition.get_el_amt_dict()
    if len(c) != 3 or sorted(c.values()) != [1, 1, 4]: return None
    x = next(e for e, n in c.items() if n == 4); a, b = [e for e in c if e != x]
    if x not in XL: return None
    for A, B in ((a, b), (b, a)):
        if A in AL and B in BL and A != B: return A, B, x
    return None

if __name__ == "__main__":
    t = time.time(); E = []; cnt = {}
    def add(src, i, s, hull):
        r = roles(s)
        if r and len(s) <= MAXN and hull is not None and hull <= HULL:
            E.append(dict(id=i, src=src, form="".join(r) + "4", hull=hull, s=s, roles=r)); cnt[src] = cnt.get(src, 0) + 1
    for f in sorted(glob.glob("discovery/cache/alex_*.json")):
        for e in json.load(open(f)):
            if "lattice_vectors" in e: add("alexandria", e["id"], struct(e), e.get("_alexandria_hull_distance"))
    seen = {e["id"] for e in E}
    for e in alex():
        if e["id"] not in seen: add("alexandria", e["id"], struct(e), e["_alexandria_hull_distance"])
    for d in mp(): add("mp", d["material_id"], Structure.from_dict(d["structure"]), d["energy_above_hull"])
    oq = "unreachable"
    try:
        r = get("https://oqmd.org/optimade/structures?" + urllib.parse.urlencode({"filter": 'chemical_formula_anonymous="A4BC" AND elements HAS ANY "F","Cl","Br","I" AND nsites<=48 AND _oqmd_stability<=0.06', "page_limit": 500}))
        for d in r["data"]:
            a = d["attributes"]; add("oqmd", "oqmd-" + str(d["id"]), Structure(Lattice(a["lattice_vectors"]), a["species_at_sites"], a["cartesian_site_positions"], coords_are_cartesian=True), a.get("_oqmd_stability"))
        oq = f"first page only ({len(r['data'])} entries returned)"
    except Exception as ex: oq = "unreachable: " + repr(ex)[:80]
    print("entries with valid roles, <=48 atoms, hull<=60:", cnt, "OQMD:", oq, f"t={time.time()-t:.0f}s", flush=True)
    for e in E:
        s = e["s"].copy(); s.replace_species({k: v for k, v in zip(e["roles"], ("Na", "Bi", "Cl")) if k != v}); e["d"] = prim(s)
        e["key"] = (len(e["d"]), SpacegroupAnalyzer(e["d"], symprec=0.1).get_space_group_number())
    old = {g[0] for g in json.load(open("discovery/cache/prototype_groups.json"))}; oldall = {i for g in json.load(open("discovery/cache/prototype_groups.json")) for i in g}
    # note: dedup only inside (n_primitive, space group at symprec 0.1) buckets; a missed duplicate just costs one extra relaxation
    reps = [e for e in E if e["id"] in old]; new = []; ndup_old = ndup_new = 0
    for e in sorted((e for e in E if e["id"] not in oldall), key=lambda e: (e["key"][0], e["hull"])):
        if any(r["key"] == e["key"] and SM.fit(r["d"], e["d"]) for r in reps): ndup_old += 1; continue
        if any(r["key"] == e["key"] and SM.fit(r["d"], e["d"]) for r in new): ndup_new += 1; continue
        new.append(e)
    out = []
    for e in new:
        s = e["d"].copy(); s.scale_lattice(28.0*len(s))
        out.append(dict(id=e["id"], src=e["src"], form=e["form"], hull=e["hull"], n=len(s), sgnum=e["key"][1], structure=s.as_dict()))
    json.dump(dict(counts=cnt, oqmd=oq, n_entries=len(E), n_old_reps_found=len(reps), n_in_old_groups=sum(e["id"] in oldall for e in E), dup_of_old=ndup_old, dup_of_new=ndup_new, protos=out), open(C/"protos.json", "w"))
    import collections
    print("old reps found", len(reps), "| duplicates of old prototypes", ndup_old, "| duplicates among new", ndup_new, "| NEW prototypes", len(out),
          dict(collections.Counter(e["src"] for e in new)), "by n:", sorted(collections.Counter(len(e["d"]) for e in new).items()), f"t={time.time()-t:.0f}s")
