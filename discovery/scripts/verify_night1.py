import re, numpy as np, spglib
from pymatgen.core import Structure, Lattice
from pymatgen.analysis.structure_matcher import StructureMatcher
RY=13605.693123
R='anish_results/ionic_discovery/dft/dft_nabicl4/'
J='overnight/night1/jobs/'
F={'ref':R+'NaBiCl4_d3.out','01':J+'01_NaBiCl4_agm088474086_vcrelax.out','02':J+'02_NaBiCl4_agm084661913_vcrelax.out'}
def last_e(f): return float(re.findall(r'^!.*=\s+(-?\d+\.\d+)',open(f).read(),re.M)[-1])
def last_h(f):
    m=re.findall(r'Final enthalpy\s+=\s+(-?\d+\.\d+)',open(f).read()); return float(m[-1]) if m else None
def last_struct(f):
    L=open(f).read().splitlines()
    ci=max(i for i,l in enumerate(L) if l.startswith('CELL_PARAMETERS'))
    ai=max(i for i,l in enumerate(L) if l.startswith('ATOMIC_POSITIONS'))
    lat=[[float(x) for x in L[ci+k].split()] for k in (1,2,3)]
    sp=[];fc=[]
    for l in L[ai+1:]:
        p=l.split()
        if len(p)<4: break
        sp.append(p[0]);fc.append([float(x) for x in p[1:4]])
    return Structure(Lattice(lat),sp,fc)
S={k:last_struct(f) for k,f in F.items()}
E={k:last_e(f) for k,f in F.items()}; H={k:last_h(f) for k,f in F.items()}
print('ctrl 00 E',last_e(J+'00_control_NaBiCl4_agm005743614_scf.out'),'ref E',E['ref'])
for k in F:
    n=len(S[k]); h=H[k]
    print(k,'E=%.8f'%E[k],'H=%.8f'%h,'E-H=%.3f mRy = %.2f meV/atom'%((E[k]-h)*1e3,(E[k]-h)*RY/n))
for k in ('01','02'):
    print('dE vs ref',k,'%.3f meV/atom'%((E[k]-E['ref'])*RY/12),' using enthalpies: %.3f'%((H[k]-H['ref'])*RY/12))
print('02-01: %.3f'%((E['02']-E['01'])*RY/12))
for k,s in S.items():
    cell=(s.lattice.matrix,s.frac_coords,[s.species[i].Z for i in range(len(s))])
    sg=[spglib.get_spacegroup(cell,symprec=p) for p in (0.01,0.1)]
    cn=lambda el:[len([n for n in s.get_neighbors(site,3.3) if n.specie.symbol=='Cl']) for site in s if site.specie.symbol==el]
    bi=[min(n.nn_distance for n in s.get_neighbors(site,4) if n.specie.symbol=='Cl') for site in s if site.specie.symbol=='Bi']
    print(k,sg,'V/atom %.3f'%(s.volume/len(s)),'V %.2f'%s.volume,'abc',np.round(s.lattice.abc,3),'angles',np.round(s.lattice.angles,2),'BiCN',cn('Bi'),'NaCN',cn('Na'),'minBiCl',np.round(bi,3))
for a,b in (('ref','01'),('ref','02'),('01','02')):
    for kw in (dict(ltol=.2,stol=.3,angle_tol=5),dict(ltol=.4,stol=.6,angle_tol=10)):
        print(a,b,kw,StructureMatcher(**kw).fit(S[a],S[b]))
# decomposition
nacl=last_e(R+'NaCl_d3.out'); bicl3=last_e(R+'BiCl3_d3.out')
nNaCl=1.0 # 2 atoms -> 1 f.u.
nBiCl3=4.0 # 16 atoms -> 4 f.u.
ref_mix=nacl/nNaCl+bicl3/nBiCl3
print('NaCl/fu',nacl/nNaCl,'BiCl3/fu',bicl3/nBiCl3)
for k in F:
    ef=E[k]/2  # 12 atoms = 2 f.u.
    d=(ef-ref_mix)*RY/6
    print(k,'decomp E (meV/atom, <0 = stable vs NaCl+BiCl3)  %.2f'%d)
# using enthalpy for NaCl/BiCl3 sensitivity
print('NaCl H-E',(-188.8728120102-nacl)*1e3,'mRy; BiCl3 H-E',(-4578.7588837712-bicl3)*1e3,'mRy')
