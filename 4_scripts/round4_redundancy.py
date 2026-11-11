import json, collections
from Bio import SeqIO
W='/home/claude/ws/audit/spdsyn/'
aln={r.id:str(r.seq) for r in SeqIO.parse('fam_trim.aln','fasta')}
orig=[r.id for r in SeqIO.parse(W+'final.fa','fasta')]
lab=json.load(open(W+'labels2.json')); cand=json.load(open('cand_all.json'))
kept=[k for k in json.load(open('screen.json'))['kept'] if not cand[k]['sp'].startswith('uncultured marine')]
def G(k): return lab[k]['grp'] if k in lab else cand[k]['grp']
def ident(a,b):
    n=m=0
    for x,y in zip(aln[a],aln[b]):
        if x!='-' and y!='-': n+=1; m+=(x==y)
    return (m/n if n>=100 else 0)
final=list(orig)
for k in sorted(kept,key=lambda k:-len(cand[k]['seq'].replace('X',''))):
    if all(not(G(o)==G(k) and ident(k,o)>=0.90) for o in final): final.append(k)
new=[k for k in final if k not in orig]
print('new kept after 90% reduction:',len(new),collections.Counter(cand[k]['grp'] for k in new))
print('species per group:',{g:len({cand[k]['sp'] for k in new if cand[k]['grp']==g}) for g in ['CRY','KAT','DISC','ALV','RHI','HAP']})
for g in ['CRY','KAT','DISC','ALV']:
    print(g, sorted({cand[k]['sp'] for k in new if cand[k]['grp']==g}))
json.dump(new,open('new_final.json','w'))
