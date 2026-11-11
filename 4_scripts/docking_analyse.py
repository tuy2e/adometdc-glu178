import sys, json, numpy as np, glob, re
from Bio.PDB import PDBParser
put=np.load('put_xyz.npy'); pc=put.mean(0)
prep=json.load(open('prep.json'))
p=PDBParser(QUIET=True)
def poses(fn):
    out=[];cur=[];sc=None
    for l in open(fn):
        if l.startswith('MODEL'): cur=[]
        elif l.startswith('REMARK VINA RESULT') or l.startswith('REMARK VINA'):
            try: sc=float(l.split()[3])
            except: pass
        elif l.startswith(('ATOM','HETATM')):
            el=l[77:79].strip() or l[12:14].strip()
            cur.append((l[12:16].strip(),el,np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
        elif l.startswith('ENDMDL'): out.append((sc,cur))
    return out
def heavy(c): return [x for x in c if x[1] in ('C','N','NA','N1') or x[0].startswith(('C','N'))]
def symrmsd(c):
    h=np.array([x[2] for x in heavy(c)])
    if len(h)!=len(put): return None
    # crystal order N1 C1 C2 C3 C4 N2 assumed linear; compare both directions after ordering by distance along chain
    a=h; best=1e9
    import itertools
    # brute force over the 2 chain directions using nearest assignment
    from scipy.optimize import linear_sum_assignment
    d=np.linalg.norm(a[:,None,:]-put[None,:,:],axis=2); r,cidx=linear_sum_assignment(d)
    return float(np.sqrt(np.mean(d[r,cidx]**2)))
def contacts(rec,c,cut=4.0):
    s=p.get_structure('r',rec)[0]; atoms=[a for a in s.get_atoms()]
    xyz=np.array([a.coord for a in atoms])
    Ns=[x[2] for x in c if x[1].startswith('N') or x[0].startswith('N')]
    hits=set()
    for n in Ns:
        d=np.linalg.norm(xyz-n,axis=1)
        for i in np.where(d<cut)[0]:
            a=atoms[i]; r=a.get_parent()
            if a.element in ('O','N'): hits.add('%s%d'%(r.get_resname(),r.id[1]))
    return sorted(hits,key=lambda x:int(re.sub(r'\D','',x)))
def cavity(rec,radius=4.0,probe=3.0,step=0.5):
    s=p.get_structure('r',rec)[0]; xyz=np.array([a.coord for a in s.get_atoms() if a.element!='H'])
    g=np.arange(-radius,radius+1e-9,step); pts=np.array([[x,y,z] for x in g for y in g for z in g]); pts=pts[np.linalg.norm(pts,axis=1)<=radius]+pc
    near=xyz[np.linalg.norm(xyz-pc,axis=1)<radius+probe+1]
    free=[np.min(np.linalg.norm(near-q,axis=1))>=probe for q in pts]
    return round(sum(free)*step**3,1)
rows=[]
for name in sys.argv[1:]:
    rec='%s_aligned.pdb'%name
    allp=[]
    for fn in sorted(glob.glob('out_%s_s*.pdbqt'%name)): allp+=poses(fn)
    allp.sort(key=lambda x:x[0])
    best=allp[0]; bc=np.array([x[2] for x in heavy(best[1])]).mean(0)
    insite=[sc for sc,c in allp if np.linalg.norm(np.array([x[2] for x in heavy(c)]).mean(0)-pc)<=2.5]
    best_in=min(insite) if insite else None
    lab=prep.get(name,{}).get('label','human 1MSV crystal' if name=='human_1MSV' else ('human E178Q (in silico)' if name=='humanE178Q' else name))
    # contacts of best in-site pose
    site_pose=[c for sc,c in allp if np.linalg.norm(np.array([x[2] for x in heavy(c)]).mean(0)-pc)<=2.5]
    cts=contacts(rec,site_pose[0]) if site_pose else []
    rows.append(dict(name=name,label=lab,best=best[0],best_dist=round(float(np.linalg.norm(bc-pc)),1),best_in_site=best_in,
                     n_in_site=len(insite),n_poses=len(allp),rmsd_cryst=round(symrmsd(site_pose[0]),2) if site_pose else None,
                     contacts=cts,cavity=cavity(rec)))
    r=rows[-1]
    print('%-26s best %5.2f (%.1f A from site) | best in human site %s | %2d/%2d poses in site | RMSD to crystal %s | cavity %5.1f A3 | contacts %s'%(
        lab[:26],r['best'],r['best_dist'],('%5.2f'%best_in) if best_in is not None else '  none',r['n_in_site'],r['n_poses'],r['rmsd_cryst'],r['cavity'],' '.join(cts)))
json.dump(rows,open('dock_summary_%s.json'%('_'.join(sys.argv[1:2]) if len(sys.argv)==2 else 'all'),'w'),indent=1)
