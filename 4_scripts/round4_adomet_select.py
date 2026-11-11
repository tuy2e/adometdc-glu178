import json, subprocess, collections
from Bio import SeqIO
R=json.load(open('/home/claude/r4/blast_hits.json'))
M=json.load(open('/home/claude/ws/audit/strameno/master_annot.json')); H=M['P17707']['seq']
cur={r.id:str(r.seq).replace('-','') for r in SeqIO.parse('/home/claude/eug/adomet/expE_add.aln','fasta')}
def pw(a,b,name):
    open(name,'w').write('>H\n%s\n>T\n%s\n'%(a,b))
    out=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000',name],capture_output=True,text=True).stdout
    d={};n=None
    for l in out.splitlines():
        if l.startswith('>'): n=l[1:]; d[n]=''
        else: d[n]+=l.strip()
    return d['H'],d['T']
def ident(a,b,name='pw_tmp.fa'):
    x,y=pw(a,b,name); n=m=0
    for p,q in zip(x,y):
        if p!='-' and q!='-': n+=1; m+=(p==q)
    return m/n if n else 0, n
groups={'Euglenida':['Euglena longa','Rhabdomonas costata','Eutreptiella gymnastica'],
        'Diplonemea':['Namystynia karyoxenos','Sulcionema specki','Hemistasia phaeocysticola','Lacrimia lanifica','Rhynchopus humris','Diplonema japonicum','Diplonema papillatum','Artemidia motanka'],
        'Centroplasthelida':['Raphidocystis contractilis','Choanocystis sp.']}
sel=[]
for sub,spp in groups.items():
    for sp in spp:
        c=[r for r in R if r['search'] in ('E','F') and r['species']==sp]
        c.sort(key=lambda r:(-(r['qcov'][1]-r['qcov'][0]), -r['best_bits']))
        best=c[0]
        seq=best['seq'].replace('*','X')
        # residues at human positions via pairwise
        h,t=pw(H,seq,'pw_%s.fa'%best['acc'])
        pos={};hp=0
        for x,y in zip(h,t):
            if x!='-': hp+=1; pos[hp]=y
        hidx=[j for j,ch in enumerate(h) if ch!='-']
        seg=t[hidx[0]:hidx[-1]+1].replace('-','')
        res=''.join(pos.get(p,'?') for p in [13,15,109,111,174,176,178,256])
        sel.append(dict(id=best['acc'],species=sp,sub=sub,search=best['search'],qcov=best['qcov'],len=len(seg),bits=best['best_bits'],evalue=best['best_evalue'],ident_human=best['best_ident'],res=res,seg=seg,n_candidates=len(c)))
# redundancy vs existing sequences of same species
for s in sel:
    if s['species']=='Eutreptiella gymnastica': s['vs_existing']=ident(cur['UPI0016DBF433'],s['seg'])
    if s['species']=='Diplonema papillatum': s['vs_existing']=ident(cur['KAJ9439767.1'],s['seg'])
print('%-28s %-11s %-16s %4s %9s %7s  res(13,15,109,111,174,176,178,256)'%('species','group','acc','len','qcov','evalue'))
for s in sel: print('%-28s %-11s %-16s %4d %9s %7s  %s %s'%(s['species'],s['sub'][:11],s['id'],s['len'],'%d-%d'%tuple(s['qcov']),s['evalue'],s['res'],('vs existing: %.2f over %d'%s['vs_existing']) if 'vs_existing' in s else ''))
json.dump(sel,open('r4_adomet_candidates.json','w'),indent=1)
