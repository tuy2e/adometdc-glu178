import json, sys, collections, itertools
from Bio import SeqIO
"""greedy redundancy reduction: within same supergroup label, drop seqs >=THR identity to a kept one (aligned columns, non-gap both)"""
aln={r.id:str(r.seq) for r in SeqIO.parse(sys.argv[1],'fasta')}
lab=json.load(open(sys.argv[2])); THR=float(sys.argv[3]); out=sys.argv[4]
must=set(sys.argv[5].split(',')) if len(sys.argv)>5 else set()
def ident(a,b):
    n=m=0
    for x,y in zip(aln[a],aln[b]):
        if x!='-' and y!='-':
            n+=1; m+= x==y
    return m/n if n else 0
ids=sorted(aln, key=lambda a:(a not in must, -len(aln[a].replace('-',''))))
kept=[]
for a in ids:
    g=lab[a]['grp']
    if a in must or all(not(lab[k]['grp']==g and ident(a,k)>=THR) for k in kept): kept.append(a)
json.dump(kept,open(out,'w')); print(len(aln),'->',len(kept))
