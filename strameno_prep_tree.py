import json,subprocess
M=json.load(open('master_annot.json')); P=json.load(open('pairmap.json')); H=M['P17707']['seq']
# re-run pairwise to extract the AdoMetDC-domain segment aligned to human 1..334
def domain(acc):
    out=open('pw/%s.fa'%acc).read()
    r=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000','pw/%s.fa'%acc],capture_output=True,text=True).stdout
    a={};n=None
    for l in r.splitlines():
        if l.startswith('>'): n=l[1:]; a[n]=''
        else: a[n]+=l.strip()
    h,t=a['H'],a['T']
    hidx=[i for i,c in enumerate(h) if c!='-']
    seg=t[hidx[0]:hidx[-1]+1].replace('-','')
    return seg
keep=[]
for a,m in M.items():
    if m['engine']!='intact': continue           # active enzymes only (dead copies analysed separately)
    if 'Plasmodium' in m['org']: continue         # giant bifunctional fusions, unreliable alignment
    if m['ident']<0.28: continue                  # too divergent to align confidently
    keep.append(a)
segs={a:domain(a) for a in keep}
keep=[a for a in keep if len(segs[a])>=250]      # drop short fragments
# redundancy filter: drop sequences >=92% identical to an already-kept one in the same group (cheap k-mer proxy via pairwise identity to human is not enough; use global identity)
def ident(s1,s2):
    import difflib
    return difflib.SequenceMatcher(None,s1,s2,autojunk=False).ratio()
final=[]
for a in sorted(keep,key=lambda x:(M[x]['group']!='Stramenopiles', x!='P17707')):
    if any(M[b]['group']==M[a]['group'] and ident(segs[a],segs[b])>0.92 for b in final): continue
    final.append(a)
print('active kept',len(keep),'-> after redundancy filter',len(final))
from collections import Counter
print(Counter(M[a]['group'] for a in final))
with open('tree_in.fa','w') as f:
    for a in final: f.write('>%s\n%s\n'%(a,segs[a]))
json.dump(final,open('tree_accs.json','w'))
