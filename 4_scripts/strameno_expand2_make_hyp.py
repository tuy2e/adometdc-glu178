import json, copy, io
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
E='/root/audit/strameno/expand/'
h=json.load(open(E+'hyp_sets.json')); d=json.load(open(E+'exp_set.json'))
t=Phylo.read(E+'expA.treefile','newick')
for cl in t.find_clades():
    if not cl.is_terminal(): cl.name=None; cl.confidence=None
allt=[x.name for x in t.get_terminals()]
OPI=[k for k,v in d.items() if v['group']=='Opisthokonta']
ST=set(h['ST'])&set(allt); DP=set(h['DP'])&set(allt); KE=set(h['KE'])&set(allt); CR=set(h['CR'])&set(allt); ALV=set(h['ALV'])&set(allt)
print('sizes ST,DP,KE,CR,ALV,OPI',len(ST),len(DP),len(KE),len(CR),len(ALV),len(OPI))
def mrca(tree,S): return tree.common_ancestor([n for n in tree.get_terminals() if n.name in S])
def is_clade(tree,S):
    c=mrca(tree,S); return set(x.name for x in c.get_terminals())==set(S)
def regraft(tree, moving, target, outg):
    tree=copy.deepcopy(tree)
    tree.root_with_outgroup([n for n in tree.get_terminals() if n.name==outg][0])
    mv=mrca(tree,moving); assert set(x.name for x in mv.get_terminals())==set(moving), 'moving not clade'
    parent=tree.get_path(mv)[-2] if len(tree.get_path(mv))>1 else tree.root
    parent.clades.remove(mv)
    if len(parent.clades)==1 and parent!=tree.root:
        # collapse unary
        gp=tree.get_path(parent)[-2] if len(tree.get_path(parent))>1 else tree.root
        only=parent.clades[0]; only.branch_length=(only.branch_length or 0)+(parent.branch_length or 0)
        gp.clades[gp.clades.index(parent)]=only
    tg=mrca(tree,target)
    tp=tree.get_path(tg)[-2] if len(tree.get_path(tg))>1 else tree.root
    new=Clade(branch_length=0.05, clades=[tg, mv])
    tp.clades[tp.clades.index(tg)]=new
    return tree
def constraint(groups, allt):
    parts=['('+','.join(sorted(g))+')' for g in groups]
    return '('+','.join(parts)+');'
out={}
GREEN=set(k for k,v in d.items() if v.get('sub') in ('land plants','Chlorophyta','other green'))&set(allt)
print('GREEN',len(GREEN))
for name,target,outg in [('li',set(OPI),sorted(GREEN)[0]),('sar',ALV,sorted(OPI)[0])]:
    tt=copy.deepcopy(t); tt.root_with_outgroup([n for n in tt.get_terminals() if n.name==outg][0]); mem=set(x.name for x in mrca(tt,target).get_terminals())
    print(name,'target clade size',len(mem),'extra members',len(mem-target),'overlap with KE|CR|GREEN|ST|DP:',len(mem&(KE|CR|GREEN|ST|DP)))
    st=regraft(t, ST|DP, target, outg)
    sio=io.StringIO(); Phylo.write(st,sio,'newick'); s=sio.getvalue().strip()
    open(f'start_{name}.nwk','w').write(s+'\n')
    open(f'h_{name}.nwk','w').write(constraint([ST|DP|target, KE|CR, GREEN], allt)+'\n')
    out[name]=True
print('ST|DP clade in ML:',is_clade(t,ST|DP),' KE|CR:',is_clade(t,KE|CR))
