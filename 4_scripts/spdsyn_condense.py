import json, sys, collections
from Bio import Phylo
lab=json.load(open('labels2.json'))
t=Phylo.read(sys.argv[1],'newick')
P=[n for n in t.get_terminals() if lab[n.name]['grp']=='PROK']
ca=t.common_ancestor(P)
if set(x.name for x in ca.get_terminals())==set(n.name for n in P): t.root_with_outgroup(ca)
else: t.root_with_outgroup(P[0])
def key(a):
    l=lab[a]; g=l['grp']; s=l['sub']
    if g=='ST': return 'ST:'+s
    if g=='ALV': return 'ALV:'+s
    if g=='DISC': return 'DISC:'+s
    return g
def walk(c, depth=0):
    tips=[x.name for x in c.get_terminals()]
    ks=collections.Counter(key(a) for a in tips)
    sup=c.name if c.name else ''
    if len(ks)==1 or len(tips)<=2:
        print('  '*depth+f"[{len(tips)}] "+', '.join(f"{k}x{v}" for k,v in ks.items())+(f"  ({sup})" if sup else '')+('  e.g. '+lab[tips[0]]['org'][:28] if len(tips)<=2 else ''))
        return
    print('  '*depth+f"+ {len(tips)} tips ({sup})")
    for ch in c.clades: walk(ch, depth+1)
walk(t.root)
