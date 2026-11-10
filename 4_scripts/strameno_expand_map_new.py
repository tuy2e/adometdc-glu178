import json,subprocess
from concurrent.futures import ThreadPoolExecutor
M=json.load(open('../master.json')); H=M['P17707']['seq']
R=json.load(open('new_recs.json'))
POS={'L13':13,'F111':111,'D174':174,'T176':176,'E178':178,'E256':256,'F285':285,'E15':15,'E67':67,'S68':68,'C82':82,'S229':229,'H243':243,'F7':7,'F223':223,'E247':247}
def run(i):
    r=R[i]; fa='pw/n%d.fa'%i
    open(fa,'w').write('>H\n%s\n>T\n%s\n'%(H,r['seq']))
    out=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000',fa],capture_output=True,text=True).stdout
    a={};n=None
    for l in out.splitlines():
        if l.startswith('>'): n=l[1:]; a[n]=''
        else: a[n]+=l.strip()
    h,t=a['H'],a['T']; res={}; hc=tc=0; first=last=None
    for x,y in zip(h,t):
        if y!='-': tc+=1
        if x!='-':
            hc+=1
            if y!='-':
                if first is None: first=tc
                last=tc
            for k,p in POS.items():
                if p==hc: res[k]=(y.upper() if y!='-' else '-', tc if y!='-' else None)
    pairs=[(x,y) for x,y in zip(h,t) if x!='-' and y!='-']
    ident=sum(1 for x,y in pairs if x.upper()==y.upper())/max(1,len(pairs))
    return i,res,round(ident,3),len(pairs),first,last
with ThreadPoolExecutor(2) as ex: out=list(ex.map(run,range(len(R))))
for i,res,ident,npair,first,last in out:
    r=R[i]; r['res']=res; r['ident']=ident; r['naligned']=npair; r['dom']=(first,last)
    g=lambda k: res.get(k,('-',None))[0]
    r['engine']='intact' if (g('S68')=='S' and g('C82')=='C' and g('H243')=='H') else ('DEAD' if (g('C82')!='C' and g('H243')!='H') else 'partial')
    r['E178']=g('E178')
json.dump(R,open('new_mapped.json','w'))
print('done',len(R))
