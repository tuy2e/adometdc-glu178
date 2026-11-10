import json, subprocess, sys
from Bio import SeqIO
M=json.load(open('/home/claude/ws/audit/strameno/master_annot.json')); H=M['P17707']['seq']
eg={r.id:str(r.seq) for r in SeqIO.parse('/home/claude/eug/euglena_gracilis_from_tblastn.fa','fasta')}
T=eg['Egr_AdoMetDC_GDJR01034331']
open('/home/claude/eug/adomet/pw_egr.fa','w').write('>H\n%s\n>T\n%s\n'%(H,T))
r=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000','/home/claude/eug/adomet/pw_egr.fa'],capture_output=True,text=True).stdout
a={};n=None
for l in r.splitlines():
    if l.startswith('>'): n=l[1:]; a[n]=''
    else: a[n]+=l.strip()
h,t=a['H'],a['T']; hidx=[j for j,c in enumerate(h) if c!='-']
seg=t[hidx[0]:hidx[-1]+1].replace('-','')
# residues at human positions
pos={}; hp=0
for j,(x,y) in enumerate(zip(h,t)):
    if x!='-':
        hp+=1; pos[hp]=(x,y)
print('domain length',len(seg),'of',len(T))
for p in [13,15,109,111,174,176,178,256]: print(p,pos[p])
open('/home/claude/eug/adomet/egr_seg.fa','w').write('>Egr_AdoMetDC_GDJR01034331\n%s\n'%seg)
