import json, sys, collections
from Bio import Phylo
sys.path.insert(0,'/home/claude/scripts')
from biparts import biparts, has
L=json.load(open('/home/claude/r4/spd/labelsR4.json'))
tf=sys.argv[1]; t=Phylo.read(tf,'newick'); allt={x.name for x in t.get_terminals()}
_,bp=biparts(tf)
G=lambda a:L[a]['grp']; SB=lambda a:L[a]['sub']
S=lambda f:{a for a in allt if f(a)}
sets={'KE':S(lambda a:SB(a)=='Kinetoplastea'),'EU':S(lambda a:SB(a)=='Euglenida'),'DIPLO':S(lambda a:SB(a)=='Diplonemea'),
 'CR':S(lambda a:G(a)=='CRY'),'KAT':S(lambda a:G(a)=='KAT'),'ST':S(lambda a:G(a)=='ST'),'DINO':S(lambda a:SB(a) in('Dinophyceae','Perkinsozoa')),
 'APIC':S(lambda a:SB(a) in('Apicomplexa','Colpodellida','Chromerida')),'HAP':S(lambda a:G(a)=='HAP'),'RHI':S(lambda a:G(a)=='RHI'),
 'RHO':S(lambda a:G(a)=='RHO'),'PROK':S(lambda a:G(a)=='PROK'),'HET':S(lambda a:SB(a)=='Heterolobosea')}
print({k:len(v) for k,v in sets.items()})
combos=[['KE'],['EU'],['DIPLO'],['KE','EU'],['EU','DIPLO'],['KE','EU','DIPLO'],['CR'],['CR','KAT'],['KE','EU','CR'],['KE','CR'],['KE','EU','DIPLO','CR'],['ST'],['ST','DINO'],['ST','DINO','KE','EU','CR'],['DINO'],['HAP'],['RHI'],['CR','RHO'],['CR','HAP']]
for c in combos:
    U=set().union(*[sets[x] for x in c]); print('%-22s %3d'%('+'.join(c),len(U)), has(allt,bp,U))
P=sets['PROK']; ca=t.common_ancestor([n for n in t.get_terminals() if n.name in P])
t.root_with_outgroup(ca if {x.name for x in ca.get_terminals()}==P else sorted(P)[0])
def smallest(names):
    ca=t.common_ancestor([n for n in t.get_terminals() if n.name in names])
    tips=[x.name for x in ca.get_terminals()]
    return len(tips), ca.confidence, dict(collections.Counter(G(a)+'/'+SB(a) for a in tips if a not in names).most_common(6))
def sister(names):
    ca=t.common_ancestor([n for n in t.get_terminals() if n.name in names]); p=t.get_path(ca); par=p[-2] if len(p)>1 else t.root
    out=collections.Counter()
    for s in [c for c in par.clades if c!=ca]:
        for x in s.get_terminals(): out[G(x.name)+'/'+SB(x.name)]+=1
    return dict(out.most_common(5)), 'parent support', par.confidence
for k in ['EU','DIPLO','KE','CR','KAT','DINO','ST','RHI','HAP']:
    print('\n',k,'| smallest clade containing all:',smallest(sets[k]),'\n    sister:',sister(sets[k]))
Phylo.write(t, tf.replace('.treefile','_rooted.nwk'),'newick')
