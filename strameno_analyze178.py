import json,csv
from collections import Counter,defaultdict
M=json.load(open('master.json')); P=json.load(open('pairmap.json'))
FT={}
for r in csv.DictReader(open('/mnt/user-data/uploads/unlabeled_fp_search/AdoMetDC_project/8_extra_runs/2_pocket/switch_calls.tsv'),delimiter='\t'):
    FT[r['candidate']]=r
def aa(acc,k): return P[acc]['res'].get(k,('-',None))[0]
# 1) concordance with FOLDTRACE structural calls
conc=Counter(); disc=[]
for acc in M:
    if acc in FT and FT[acc]['fold_ok']=='true':
        for k,ftk in [('E178','Glu178'),('D174','Asp174'),('E256','Glu256'),('L13','Leu13'),('F111','Phe111'),('F285','Phe285'),('T176','Thr176')]:
            f=FT[acc][ftk+'_obs']; f=f[0] if f not in ('-','') else '-'
            s=aa(acc,k)
            if f=='-' or s=='-': conc['n/a']+=1; continue
            if f==s: conc['agree']+=1
            else: conc['disagree']+=1; disc.append((acc,M[acc]['org'][:25],k,'FT='+FT[acc][ftk+'_obs'],'seq='+s+str(P[acc]['res'][k][1])))
print('FOLDTRACE-structure vs sequence-alignment agreement (pocket positions):',dict(conc))
for d in disc: print('   disagree:',d)
# 2) engine status
for acc,m in M.items():
    c82,h243,s68=aa(acc,'C82'),aa(acc,'H243'),aa(acc,'S68')
    m['engine']='intact' if (s68=='S' and c82=='C' and h243=='H') else ('DEAD' if (c82!='C' and h243!='H') else 'partial')
    m['E178']=aa(acc,'E178'); m['ident']=P[acc]['ident']
json.dump(M,open('master_annot.json','w'))
print('\nEngine status:',Counter(m['engine'] for m in M.values()))
print('dead copies:',[(a,m['org'][:28],m['group']) for a,m in M.items() if m['engine']=='DEAD'])
print('partial engine:',[(a,m['org'][:28],m['group'],'S68='+aa(a,'S68'),'C82='+aa(a,'C82'),'H243='+aa(a,'H243')) for a,m in M.items() if m['engine']=='partial'])
# 3) residue at 178 by group (active enzymes only)
print('\nResidue at human-178 position, ACTIVE (intact-engine) enzymes:')
tab=defaultdict(Counter)
for m in M.values():
    if m['engine']!='intact': continue
    tab[(m['group'],m['sub'])][m['E178']]+=1
for k in sorted(tab): print('  %-15s %-20s %s'%(k[0],k[1],dict(tab[k])))
