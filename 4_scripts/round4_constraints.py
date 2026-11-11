import json, re
from Bio import SeqIO
# --- SpdSyn hypotheses
L=json.load(open('/home/claude/r4/spd/labelsR4.json'))
ids=[r.id for r in SeqIO.parse('/home/claude/r4/spd/spdR4.trim.aln','fasta')]
G=lambda k:L[k]['grp']; SB=lambda k:L[k]['sub']
S=lambda f:[k for k in ids if f(k)]
sets={'KE':S(lambda k:SB(k)=='Kinetoplastea'),'EU':S(lambda k:SB(k)=='Euglenida'),'DIPLO':S(lambda k:SB(k)=='Diplonemea'),
 'CR':S(lambda k:G(k)=='CRY'),'ST':S(lambda k:G(k)=='ST'),'DINO':S(lambda k:SB(k) in ('Dinophyceae','Perkinsozoa')),'RHO':S(lambda k:G(k)=='RHO')}
print({k:len(v) for k,v in sets.items()})
H={'h1_crkeeu':['CR','KE','EU'],'h1b_crke':['CR','KE'],'h1c_creuglenozoa':['CR','KE','EU','DIPLO'],'h2_crrho':['CR','RHO'],'h3_stdino':['ST','DINO'],'h4_lineage':['ST','DINO','KE','EU','CR']}
for h,gs in H.items():
    g=sorted(set().union(*[sets[x] for x in gs])); rest=[k for k in ids if k not in g]
    open('/home/claude/r4/spd/%s.nwk'%h,'w').write('(('+','.join(g)+'),'+','.join(rest)+');\n')
    print(h,len(g),'constrained')
json.dump(sets,open('/home/claude/r4/spd/sets.json','w'))
# --- AdoMetDC H_kin with the new euglenozoans
s=open('/home/claude/ws/audit/strameno/expand/h_kin.nwk').read().strip()
grs=re.findall(r'\(([^()]+)\)',s)
top=s[1:-2]
for g in grs: top=top.replace('('+g+')','',1)
rest=[x for x in top.split(',') if x]
new=json.load(open('/home/claude/r4/adomet/new12_labels.json'))
eugleno_new=['Egr_AdoMetDC_GDJR01034331']+[k for k,v in new.items() if v['sub'] in ('Euglenida','Diplonemea')]
centro=[k for k,v in new.items() if v['sub']=='Centrohelida']
g68=[x for x in grs if 'UPI0016DBF433' not in x][0].split(','); g8=[x for x in grs if 'UPI0016DBF433' in x][0].split(',')
ids_a=[r.id for r in SeqIO.parse('/home/claude/r4/adomet/r4_trim.aln','fasta')]
out='(('+','.join(g68)+'),('+','.join(g8+eugleno_new)+'),'+','.join(rest+centro)+');\n'
allt=set(re.findall(r'[A-Za-z0-9_.]+',out)); assert allt==set(ids_a), (len(allt),len(ids_a),set(ids_a)-allt)
open('/home/claude/r4/adomet/h_kin_R4.nwk','w').write(out)
print('H_kin R4: Euglenozoa group',len(g8+eugleno_new),'| rest-of-lineage group',len(g68),'| free',len(rest+centro))
