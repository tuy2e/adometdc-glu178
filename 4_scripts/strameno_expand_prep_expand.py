import json, subprocess, numpy as np, collections
M=json.load(open('../master_annot.json')); H=M['P17707']['seq']
R={r['id']:(i,r) for i,r in enumerate(json.load(open('new_mapped.json')))}
S=json.load(open('screen_recs.json')); ident=np.load('ident.npy'); sidx={r['id']:i for i,r in enumerate(S)}
old_tree=json.load(open('../tree_accs.json'))
# old domain segments
oldseg={}; k=None
for l in open('../tree_in.fa'):
    l=l.strip()
    if l.startswith('>'): k=l[1:]; oldseg[k]=''
    else: oldseg[k]+=l
def subname(r):
    lin=r['lineage']; q=r['q']
    if r['sup']=='Stramenopiles': return 'Bigyra' if 'Bigyra' in lin else ('Ochrophyta' if 'Ochrophyta' in lin else ('Oomycota' if 'Oomycota' in lin else 'other stramenopile'))
    if q in ('StramOther',): return 'other stramenopile'
    if q=='KinetoOther': return 'Kinetoplastea'
    if q=='Rhodophyta': return 'Rhodophyta'
    if q=='Glaucophyta': return 'Glaucophyta'
    if q=='Apusozoa': return 'Apusozoa'
    if r['sup']=='Alveolata' and q=='Dinophyceae':
        return 'Dinophyceae (ODC fusion)' if len(r['seq'])>=700 else 'Dinophyceae'
    return q
def domain(i):
    r=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000','pw/n%d.fa'%i],capture_output=True,text=True).stdout
    a={};n=None
    for l in r.splitlines():
        if l.startswith('>'): n=l[1:]; a[n]=''
        else: a[n]+=l.strip()
    h,t=a['H'],a['T']; hidx=[j for j,c in enumerate(h) if c!='-']
    return t[hidx[0]:hidx[-1]+1].replace('-','')
cands={}; dropped=[]
for rid,(i,r) in R.items():
    if r['engine']!='intact' or r['ident']<0.28: continue
    if r['sup'].startswith('MULTI') or r['sup']=='BACTERIA': dropped.append((rid,r['org'][:40],'multi-organism record')); continue
    si=sidx[rid]; row=ident[si]
    foreign=[j for j in range(len(S)) if S[j]['sup']!=r['sup']]
    own=[j for j in range(len(S)) if S[j]['sup']==r['sup'] and j!=si]
    bf=max(row[foreign]) if foreign else 0; bo=max(row[own]) if own else 0
    if bf>=0.50 and bf>bo: dropped.append((rid,r['org'][:40],'closest relative in another supergroup (%.2f vs own %.2f)'%(bf,bo))); continue
    seg=domain(i)
    if len(seg)<250: dropped.append((rid,r['org'][:40],'domain <250 aa')); continue
    org=r['org'].split(';')[0]
    cands[rid]=dict(id=rid,org=org,group=r['sup'] if not r['sup'].startswith('other') else 'Apusozoa',sub=subname(r),seg=seg,src=r['src'],E178=r['E178'],new=True)
for a in old_tree:
    cands[a]=dict(id=a,org=M[a]['org'],group=M[a]['group'],sub=M[a]['sub'],seg=oldseg[a],src='UniProt/AF',E178=M[a]['E178'],new=False)
print('candidates',len(cands),'dropped',len(dropped))
for d in dropped: print('  DROP',d)
# diversity-capped selection per subgroup
CAP=collections.defaultdict(lambda:99)
CAP.update({'Metazoa':14,'Fungi':10,'land plants':9,'Chlorophyta':10,'Ochrophyta':22,'Oomycota':12,'Haptophyta':12,'Rhodophyta':12,'Rhizaria':12,'Dinophyceae':16,'Dinophyceae (ODC fusion)':8,'Heterolobosea':8,'other opisthokonts':5})
REF={'P17707','P21182','Q04694','O76240','Q587A7','O02655','B7FPB8','M4BIM5','D7FY29','A0A5A8CY83','L1JMA4','UPI0016DBF433','L8GRS6','Q96286','C5LKE7'}
def I(a,b):
    if a in sidx and b in sidx: return ident[sidx[a],sidx[b]]
    import difflib; return difflib.SequenceMatcher(None,cands[a]['seg'],cands[b]['seg'],autojunk=False).ratio()
bysub=collections.defaultdict(list)
for a,c in cands.items(): bysub[(c['group'],c['sub'])].append(a)
final=[]
for key,mem in sorted(bysub.items()):
    sel=[a for a in mem if a in REF]
    pool=[a for a in mem if a not in sel]
    if not sel: sel=[max(pool,key=lambda a:len(cands[a]['seg']))]; pool.remove(sel[0])
    while pool and len(sel)<CAP[key[1]]:
        nxt=min(pool,key=lambda a:max(I(a,b) for b in sel))
        if max(I(nxt,b) for b in sel)>0.92: break
        sel.append(nxt); pool.remove(nxt)
    final+=sel
    print('%-16s %-26s %3d -> %3d   178: %s'%(key[0][:16],key[1][:26],len(mem),len(sel),dict(collections.Counter(cands[a]['E178'] for a in sel))))
print('TOTAL',len(final))
json.dump({a:{k:v for k,v in cands[a].items()} for a in final},open('exp_set.json','w'))
with open('exp_in.fa','w') as f:
    for a in final: f.write('>%s\n%s\n'%(a,cands[a]['seg']))
json.dump(dropped,open('exp_dropped.json','w'))
