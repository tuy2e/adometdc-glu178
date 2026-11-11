import sys, json, numpy as np, re
sys.path.insert(0,'/home/claude/scripts')
import importlib.util
spec=importlib.util.spec_from_file_location('da','/home/claude/scripts/dock_analyse.py')
src=open('/home/claude/scripts/dock_analyse.py').read().split('rows=[]')[0]
exec(src)
names=sys.argv[1:]; rows=[]
for name in names:
    rec='%s_aligned.pdb'%name
    lab=prep.get(name,{}).get('label','human 1MSV crystal' if name=='human_1MSV' else ('human E178Q (in silico)' if name=='humanE178Q' else name))
    def summ(fn):
        P=poses(fn); P.sort(key=lambda x:x[0])
        ins=[(sc,c) for sc,c in P if np.linalg.norm(np.array([x[2] for x in heavy(c)]).mean(0)-pc)<=2.5]
        best=P[0]; bd=float(np.linalg.norm(np.array([x[2] for x in heavy(best[1])]).mean(0)-pc))
        return best[0],round(bd,1),(ins[0][0] if ins else None),len(ins),len(P),(contacts(rec,ins[0][1]) if ins else contacts(rec,best[1]))
    v=summ('site_%s.pdbqt'%name); a=summ('ad4out_%s.pdbqt'%name)
    rows.append(dict(name=name,label=lab,vina_site_best=v[0],vina_site_dist=v[1],vina_site_n_in=v[3],vina_contacts=v[5],
                     ad4_best=a[0],ad4_best_dist=a[1],ad4_best_in_site=a[2],ad4_n_in=a[3],ad4_contacts=a[5]))
    print('%-26s | 12A-box Vina: best %5.2f at %4.1f A, %d/9 in site | AD4 (charge): best %6.2f at %4.1f A, in-site best %s, %d/9 in site | AD4 contacts: %s'%(
        lab[:26],v[0],v[1],v[3],a[0],a[1],('%6.2f'%a[2]) if a[2] is not None else '  none',a[3],' '.join(a[5])))
json.dump(rows,open('dock_summary_site_ad4.json','w'),indent=1)
