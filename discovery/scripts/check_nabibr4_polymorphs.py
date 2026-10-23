"""Is our 'NaBiBr4 on the NaBiCl4 prototype' model already an Alexandria entry? Also list NaBiCl4 polymorphs."""
import json, sys, urllib.request, urllib.parse, glob
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.structure_matcher import StructureMatcher
U='https://alexandria.icams.rub.de/pbe/v1/structures'
def fetch(els, formula):
    f='elements HAS ALL '+','.join(f'"{e}"' for e in els)+f' AND nelements={len(els)} AND chemical_formula_reduced="{formula}"'
    r=json.load(urllib.request.urlopen(U+'?'+urllib.parse.urlencode({'filter':f,'page_limit':50}),timeout=60))
    out=[]
    for d in r['data']:
        a=d['attributes']
        s=Structure(Lattice(a['lattice_vectors']),a['species_at_sites'],a['cartesian_site_positions'],coords_are_cartesian=True)
        out.append((d['id'],a.get('_alexandria_hull_distance'),a.get('_alexandria_space_group'),s))
    return out
sm=StructureMatcher(ltol=0.3,stol=0.5,angle_tol=10)
def load(pat):
    p=glob.glob(pat); return Structure.from_file(p[0]) if p else None
mine={'Br_on_Cl_proto':load('discovery/results/mace/REF-NaBiBr4-on-NaBiCl4proto/*.cif'),'Alex_NaBiBr4_input':load('discovery/results/mace/NBCB-4-NaBiBr4/*.cif'),
      'Cl_parent':load('discovery/results/mace/NBCB-0-parent/*.cif')}
print({k:(v.composition.reduced_formula,len(v)) if v else None for k,v in mine.items()})
for els,form in [(('Na','Bi','Br'),'BiBr4Na'),(('Na','Bi','Cl'),'BiCl4Na')]:
    for id_,h,sg,s in sorted(fetch(els,form),key=lambda x:x[1] if x[1] is not None else 9):
        m={k:(sm.fit(v,s) if v.composition.reduced_formula==s.composition.reduced_formula else sm.fit_anonymous(v,s)) for k,v in mine.items() if v}
        print(form,id_,'hull_meV',None if h is None else round(h*1000,1),'sg',sg,'nat',len(s),m)
