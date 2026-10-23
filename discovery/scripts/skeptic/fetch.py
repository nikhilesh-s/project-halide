"""Fetch Alexandria PBE entries (OPTIMADE) for the skeptic checks; cache to discovery/results/skeptic/alex.json."""
import json, urllib.request, urllib.parse
U='https://alexandria.icams.rub.de/pbe/v1/structures'
FORMS=['BiBr4Na','BiCl4Na','Br4NaSb','Cl4NaSb','BiI4Na','BrNa','BiBr3','I4NaSb','Br4KBi','BiBr4Li','BiBr4K','AgBiBr4','BiBr4Rb','BiBr6Na3']
out={}
for f in FORMS:
    q=urllib.parse.urlencode({'filter':f'chemical_formula_reduced="{f}"','page_limit':100})
    r=json.load(urllib.request.urlopen(U+'?'+q,timeout=90))
    out[f]=[{'id':d['id'],**{k:v for k,v in d['attributes'].items() if k in('lattice_vectors','species_at_sites','cartesian_site_positions','nsites') or k.startswith('_alexandria')}} for d in r['data']]
    print(f,len(out[f]),[(e['id'],e.get('_alexandria_space_group'),round(1000*e['_alexandria_hull_distance'],1),e['nsites']) for e in out[f]])
json.dump(out,open('discovery/results/skeptic/alex.json','w'))
print(sorted(k for k in out['BiBr4Na'][0] if k.startswith('_alex')))
