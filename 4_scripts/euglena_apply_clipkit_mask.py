# apply an original ClipKIT column mask (from its -l log) to an alignment with the same column count
import sys
from Bio import SeqIO
log,aln,out=sys.argv[1:4]
keep=[l.split()[1]=='keep' for l in open(log) if l.strip()]
recs=list(SeqIO.parse(aln,'fasta')); L=len(recs[0].seq); assert L==len(keep),(L,len(keep))
with open(out,'w') as f:
    for r in recs: f.write('>%s\n%s\n'%(r.id,''.join(c for c,k in zip(str(r.seq),keep) if k)))
print(out,len(recs),'seqs',sum(keep),'cols')
