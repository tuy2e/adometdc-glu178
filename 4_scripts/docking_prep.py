"""Superimpose each model on 1MSV chain A around the putrescine site (CA atoms of human residues within 18 A of PUT,
mapped by a pairwise MAFFT alignment), write the receptor in the 1MSV frame, and report residues at key positions."""
import warnings; warnings.filterwarnings('ignore')
import json, subprocess, numpy as np
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
U='/mnt/user-data/uploads/unlabeled_fp_search/unlabeled-fp-search/aim1_workflow/runs/stuy__adometdc/'
p=PDBParser(QUIET=True)
H=p.get_structure('h',U+'1MSV.pdb')[0]
chA=H['A']
put=[r for r in chA if r.get_resname()=='PUT']
if not put:  # PUT may be a separate chain record
    put=[r for ch in H for r in ch if r.get_resname()=='PUT']
put=put[0]; putxyz=np.array([a.coord for a in put])
three={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','PYR':'S'}
M=json.load(open('/home/claude/ws/audit/strameno/master_annot.json')); HS=M['P17707']['seq']
# human residues of chain A (+ chain B if the alpha chain is separate)
hres={}
for ch in H:
    for r in ch:
        if r.id[0]==' ' and 'CA' in r: hres.setdefault(r.id[1],r)
def align_map(tseq):
    open('tmp.fa','w').write('>H\n%s\n>T\n%s\n'%(HS,tseq))
    out=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000','tmp.fa'],capture_output=True,text=True).stdout
    a={};n=None
    for l in out.splitlines():
        if l.startswith('>'): n=l[1:]; a[n]=''
        else: a[n]+=l.strip()
    m={};hc=tc=0
    for x,y in zip(a['H'],a['T']):
        if x!='-': hc+=1
        if y!='-': tc+=1
        if x!='-' and y!='-': m[hc]=tc
    return m
POS=[13,15,80,109,111,113,174,176,178,256,285]
models={'B7FPB8':'diatom Phaeodactylum','D7FY29':'brown alga Ectocarpus','A0A5A8CY83':'Cafeteria (Bigyra)','M4BIM5':'oomycete Hyaloperonospora',
 'T0REK4':'oomycete Saprolegnia','A0A024TG05':'oomycete Aphanomyces','L1JMA4':'cryptophyte Guillardia','C5LKE7':'Perkinsus (Asp copy)',
 'A0A2V2WHV6':'T. cruzi','L8GRS6':'Acanthamoeba','A0AAJ7UEA7':'lamprey Petromyzon','D8S2S4':'plant Selaginella','A0A7S3TQM9':'haptophyte Emiliania'}
out={}
for acc,label in models.items():
    s=p.get_structure(acc,U+'models/%s.pdb'%acc)[0]; ch=list(s)[0]
    tres=[r for r in ch if r.id[0]==' ']
    tseq=''.join(three.get(r.get_resname(),'X') for r in tres)
    m=align_map(tseq); tby={i+1:r for i,r in enumerate(tres)}
    near=[n for n,r in hres.items() if np.min(np.linalg.norm(putxyz-r['CA'].coord,axis=1))<18 and n in m and m[n] in tby and 'CA' in tby[m[n]]]
    sup=Superimposer(); sup.set_atoms([hres[n]['CA'] for n in near],[tby[m[n]]['CA'] for n in near])
    # iterative core: drop residues >2 A and refit (twice)
    for _ in range(2):
        rot,tran=sup.rotran
        d=[np.linalg.norm(hres[n]['CA'].coord-(np.dot(tby[m[n]]['CA'].coord,rot)+tran)) for n in near]
        core=[n for n,x in zip(near,d) if x<2.0]
        if len(core)>=20: sup=Superimposer(); sup.set_atoms([hres[n]['CA'] for n in core],[tby[m[n]]['CA'] for n in core]); near_core=core
    rot,tran=sup.rotran
    for atom in s.get_atoms(): atom.transform(rot,tran)
    io=PDBIO(); io.set_structure(s); io.save('%s_aligned.pdb'%acc)
    res={p_:('%s%d'%(three.get(tby[m[p_]].get_resname(),'X'),tby[m[p_]].id[1]) if p_ in m and m[p_] in tby else '-') for p_ in POS}
    pl=[a.bfactor for a in s.get_atoms() if a.get_id()=='CA']
    site_pl=np.mean([tby[m[p_]]['CA'].bfactor for p_ in POS if p_ in m and m[p_] in tby])
    out[acc]=dict(label=label,core_ca=len(near_core),rmsd=round(float(sup.rms),2),res=res,site_plddt=round(float(site_pl),1))
    print('%-11s %-27s core %3d CA, RMSD %.2f, site pLDDT %.0f | %s'%(acc,label,len(near_core),sup.rms,site_pl,' '.join('%d:%s'%(k,v) for k,v in res.items())))
json.dump(out,open('prep.json','w'),indent=1)
np.save('put_xyz.npy',putxyz)
print('PUT atoms',len(putxyz),'centroid',putxyz.mean(0).round(2))
