import json, collections, sys
from Bio import Phylo
W='/home/claude/ws/audit/spdsyn/'
t=Phylo.read('fam.ft.nwk','newick')
parent={}; depth={}
def walk(c,d):
    depth[c]=d
    for ch in c.clades: parent[ch]=c; walk(ch,d+(ch.branch_length or 0))
sys.setrecursionlimit(10000); walk(t.root,0)
tips={x.name:x for x in t.get_terminals()}
clade=set(json.load(open(W+'c1_spdsyn_clade.json')))
U={u['acc']:u for u in json.load(open(W+'uniprot_filtered.json'))}
refs=[n for n in tips if n in U or n in clade or n in json.load(open(W+'ncbi_named_keep.json'))]
cand=json.load(open('cand_all.json'))
def anc(n):
    s=[]; 
    while n in parent: n=parent[n]; s.append(n)
    return s
out={}
for cid in cand:
    if cid not in tips: continue
    a=tips[cid]; A=anc(a); Aset=set(A)
    ds=[]
    for r in refs:
        b=tips[r]; x=b
        while x not in Aset and x in parent: x=parent[x]
        ds.append((depth[a]+depth[b]-2*depth[x], r))
    ds.sort(); nn=ds[:3]
    lab=lambda r: 'SPD' if r in clade else (U[r]['name'][:30] if r in U else 'ncbi_other')
    labs=[lab(r) for _,r in nn]
    isspd = labs[0]=='SPD' and labs.count('SPD')>=2
    out[cid]=dict(spd=isspd, nn=[(round(d,3),r,l) for (d,r),l in zip(nn,labs)])
json.dump(out,open('classify.json','w'),indent=0)
c=collections.Counter((cand[k]['grp'],v['spd']) for k,v in out.items())
print(sorted(c.items()))
print('non-SpdSyn nearest-label counts:',collections.Counter(v['nn'][0][2] for v in out.values() if not v['spd']).most_common(8))
