import json, collections, sys
from Bio import Phylo
W='/home/claude/ws/audit/spdsyn/'
U=json.load(open(W+'uniprot_filtered.json'))
lab={u['acc']:dict(g=u['g'],org=u['org'],name=u['name'],lin=u['lin']) for u in U}
nk=json.load(open(W+'ncbi_named_keep.json'))
for k,v in nk.items(): lab[k]=dict(g=v.get('g','?') if isinstance(v,dict) else '?',org=(v.get('org','?') if isinstance(v,dict) else str(v)),name=(v.get('name','?') if isinstance(v,dict) else '?'),lin='')
clade=set(json.load(open(W+'c1_spdsyn_clade.json')))
t=Phylo.read(sys.argv[1],'newick')
tips={x.name:x for x in t.get_terminals()}
for nm in tips:
    if nm not in lab: lab[nm]=dict(g='EGR' if nm.startswith('Egr') else '?',org='Euglena gracilis' if nm.startswith('Egr') else '?',name=nm,lin='')
# root on human spermine synthase P52788? no: root at midpoint for neighbour listing only
egr=[n for n in tips if n.startswith('Egr')]
def biparts(tree):
    allt=frozenset(tips); out=[]
    for c in tree.find_clades():
        s=frozenset(x.name for x in c.get_terminals())
        if 0<len(s)<len(allt): out.append((s,c.confidence))
    return allt,out
allt,bp=biparts(t)
def is_clade(S):
    S=frozenset(S)
    for s,conf in bp:
        if s==S or (allt-s)==S: return True,conf
    return False,None
print('SpdSyn clade (c1 def, %d) monophyletic in new tree:'%len(clade&set(tips)), is_clade(clade&set(tips)))
for e in egr:
    print('\n==',e,' + SpdSyn clade monophyletic:', is_clade((clade&set(tips))|{e}))
    d=sorted(((t.distance(tips[e],tips[o]),o) for o in tips if o!=e))[:8]
    for dist,o in d: print('   %.3f %-14s %-4s %-32s %s'%(dist,o,lab[o]['g'],lab[o]['org'][:32],lab[o]['name'][:45]))
    # smallest clade containing e with >=4 tips: composition of names
    best=None
    for s,conf in bp:
        for S in (s,allt-s):
            if e in S and 5<=len(S)<=60 and (best is None or len(S)<len(best[0])): best=(S,conf)
    S,conf=best
    print('   smallest 5+ tip group (%d tips, support %s):'%(len(S),conf), dict(collections.Counter(lab[o]['name'][:30] for o in S if o!=e).most_common(5)))
