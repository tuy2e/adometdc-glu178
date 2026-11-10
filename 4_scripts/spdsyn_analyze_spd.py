import json, sys, collections
from Bio import Phylo
sys.path.insert(0,'/root/audit/strameno/expand2')
from biparts import biparts, has
lab=json.load(open('labels2.json'))
tf=sys.argv[1]
t=Phylo.read(tf,'newick'); allt=set(x.name for x in t.get_terminals())
_,bp=biparts(tf)
G=lambda a: lab[a]['grp']; SB=lambda a: lab[a]['sub']
def S(f): return {a for a in allt if f(a)}
sets={
 'ST':S(lambda a:G(a)=='ST'),
 'ST_noBigyra':S(lambda a:G(a)=='ST' and SB(a)!='Bigyra'),
 'KE':S(lambda a:SB(a)=='Kinetoplastea'),
 'CR':S(lambda a:G(a)=='CRY'),
 'DIPLO':S(lambda a:SB(a)=='Diplonemea'),
 'HET':S(lambda a:SB(a)=='Heterolobosea'),
 'DINO':S(lambda a:SB(a) in ('Dinophyceae','Perkinsozoa')),
 'APIC':S(lambda a:SB(a) in ('Apicomplexa','Colpodellida','Chromerida')),
 'HAP':S(lambda a:G(a)=='HAP'),'RHO':S(lambda a:G(a)=='RHO'),'VIR':S(lambda a:G(a)=='VIR'),
 'OPI':S(lambda a:G(a)=='OPI'),'AMO':S(lambda a:G(a)=='AMO'),'PROK':S(lambda a:G(a)=='PROK'),'RHI':S(lambda a:G(a)=='RHI'),
}
print({k:len(v) for k,v in sets.items()})
combos=[['KE','CR'],['KE','DIPLO'],['KE','CR','DIPLO'],['KE','DIPLO','HET'],['ST'],['ST_noBigyra'],['DINO'],['ST','DINO'],['ST','RHO'],['ST','KE','CR'],['ST','DINO','KE','CR'],['CR','HAP'],['DINO','APIC'],['KE'],['CR']]
for c in combos:
    U=set().union(*[sets[x] for x in c])
    print('+'.join(c), len(U), has(allt,bp,U))
# root on prokaryotes and report sister of CR, KE, ST clades
pk=sorted(sets['PROK']); t.root_with_outgroup(t.common_ancestor([n for n in t.get_terminals() if n.name in sets['PROK']]) if set(x.name for x in t.common_ancestor([n for n in t.get_terminals() if n.name in sets['PROK']]).get_terminals())==sets['PROK'] else pk[0])
def sister(names):
    ca=t.common_ancestor([n for n in t.get_terminals() if n.name in names]); p=t.get_path(ca); par=p[-2] if len(p)>1 else t.root
    sib=[c for c in par.clades if c!=ca]
    out=collections.Counter(); 
    for s in sib:
        for x in s.get_terminals(): out[G(x.name)+'/'+SB(x.name)]+=1
    return len(ca.get_terminals()), dict(out.most_common(6)), par.confidence if hasattr(par,'confidence') else None, ca.name
for k in ['CR','KE','ST','DINO','HAP']:
    if sets[k]: print('sister of',k, sister(sets[k]))
