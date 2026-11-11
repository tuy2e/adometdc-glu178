import json, collections
from Bio import SeqIO
W='/home/claude/ws/audit/spdsyn/'
aln={r.id:str(r.seq) for r in SeqIO.parse('fam_trim.aln','fasta')}
orig=[r.id for r in SeqIO.parse(W+'final.fa','fasta')]
lab=json.load(open(W+'labels2.json'))
cand=json.load(open('cand_all.json')); cl=json.load(open('classify.json'))
new=[k for k,v in cl.items() if v['spd']]
print('original in family aln:',sum(o in aln for o in orig),'/',len(orig))
def G(k): return lab[k]['grp'] if k in lab else cand[k]['grp']
SUP={'KAT':'CRY'}
def S(k): g=G(k); return SUP.get(g,g)
def ident(a,b):
    n=m=0
    for x,y in zip(aln[a],aln[b]):
        if x!='-' and y!='-': n+=1; m+=(x==y)
    return (m/n if n>=100 else 0)
pool=[o for o in orig if o in aln]+new
drop={}
for k in new:
    best_f=max(((ident(k,o),o) for o in pool if o!=k and S(o)!=S(k)),default=(0,None))
    best_o=max(((ident(k,o),o) for o in pool if o!=k and S(o)==S(k)),default=(0,None))
    sp=cand[k]['sp']
    if sp.split()[0] in ('Kryptoperidinium','Durinskia'): drop[k]='dinotom (diatom-derived endosymbiont genes)'
    elif best_f[0]>=0.50 and best_f[0]>best_o[0]:
        drop[k]='closest relative in another supergroup: %.2f %s (%s) vs own %.2f'%(best_f[0],best_f[1],G(best_f[1]),best_o[0])
for k,v in drop.items(): print('DROP',k,cand[k]['sp'],'|',v)
kept=[k for k in new if k not in drop]
print('kept',len(kept),collections.Counter(cand[k]['grp'] for k in kept))
json.dump(dict(kept=kept,drop=drop),open('screen.json','w'),indent=1)
