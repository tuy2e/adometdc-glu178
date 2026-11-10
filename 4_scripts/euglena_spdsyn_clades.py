"""Clade checks on the SpdSyn tree with Euglena (adapted from audit/spdsyn/analyze_spd.py)."""
import json, sys, collections
from Bio import Phylo
sys.path.insert(0,'/home/claude/scripts')
from biparts import biparts, has
lab=json.load(open('/home/claude/eug/spd/labels2E.json'))
tf=sys.argv[1]
t=Phylo.read(tf,'newick'); allt=set(x.name for x in t.get_terminals())
_,bp=biparts(tf)
G=lambda a: lab[a]['grp']; SB=lambda a: lab[a]['sub']
def S(f): return {a for a in allt if f(a)}
sets={'ST':S(lambda a:G(a)=='ST'),'KE':S(lambda a:SB(a)=='Kinetoplastea'),'EU':S(lambda a:SB(a)=='Euglenida'),
 'CR':S(lambda a:G(a)=='CRY'),'DIPLO':S(lambda a:SB(a)=='Diplonemea'),'HET':S(lambda a:SB(a)=='Heterolobosea'),
 'DINO':S(lambda a:SB(a) in ('Dinophyceae','Perkinsozoa')),'APIC':S(lambda a:SB(a) in ('Apicomplexa','Colpodellida','Chromerida')),
 'HAP':S(lambda a:G(a)=='HAP'),'RHO':S(lambda a:G(a)=='RHO'),'PROK':S(lambda a:G(a)=='PROK')}
print({k:len(v) for k,v in sets.items()})
combos=[['KE'],['KE','EU'],['EU','DIPLO'],['KE','EU','DIPLO'],['KE','DIPLO'],['KE','EU','CR'],['KE','CR'],['KE','EU','DIPLO','CR'],['CR'],['ST','DINO','KE','EU','CR'],['ST','DINO']]
for c in combos:
    U=set().union(*[sets[x] for x in c]); print('%-22s %3d'%('+'.join(c),len(U)), has(allt,bp,U))
P=sets['PROK']; ca=t.common_ancestor([n for n in t.get_terminals() if n.name in P])
t.root_with_outgroup(ca if set(x.name for x in ca.get_terminals())==P else sorted(P)[0])
def sister(names):
    ca=t.common_ancestor([n for n in t.get_terminals() if n.name in names]); p=t.get_path(ca); par=p[-2] if len(p)>1 else t.root
    out=collections.Counter()
    for s in [c for c in par.clades if c!=ca]:
        for x in s.get_terminals(): out[G(x.name)+'/'+SB(x.name)]+=1
    return len(ca.get_terminals()), ca.confidence if ca.confidence is not None else ca.name, dict(out.most_common(6)), 'parent support', par.confidence if par.confidence is not None else par.name
for k in ['EU','KE','CR','DIPLO']:
    print('sister of',k,sister(sets[k]))
# walk up from Euglena: first 4 enclosing clades
e=[x for x in t.get_terminals() if SB(x.name)=='Euglenida'][0]
for c in list(reversed(t.get_path(e)[:-1]))[:4]:
    tips=[x.name for x in c.get_terminals()]
    print('  enclosing clade %3d tips, support %s:'%(len(tips),c.confidence if c.confidence is not None else c.name), dict(collections.Counter(G(a)+'/'+SB(a) for a in tips).most_common(6)))
Phylo.write(t, tf.replace('.treefile','_rooted.nwk'),'newick')
