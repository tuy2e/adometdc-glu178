import json,subprocess,sys,os
from concurrent.futures import ThreadPoolExecutor
M=json.load(open('master.json')); H=M['P17707']['seq']
POS={'L13':13,'F111':111,'D174':174,'T176':176,'E178':178,'E256':256,'F285':285,'E15':15,'E67':67,'S68':68,'C82':82,'S229':229,'H243':243,'F7':7,'F223':223,'E247':247}
def run(acc):
    fa='pw/%s.fa'%acc
    open(fa,'w').write('>H\n%s\n>T\n%s\n'%(H,M[acc]['seq']))
    out=subprocess.run(['mafft','--quiet','--localpair','--maxiterate','1000',fa],capture_output=True,text=True).stdout
    a={};n=None
    for l in out.splitlines():
        if l.startswith('>'): n=l[1:]; a[n]=''
        else: a[n]+=l.strip()
    h,t=a['H'],a['T']; res={}
    hc=0; tc=0
    for i,(x,y) in enumerate(zip(h,t)):
        if y!='-': tc+=1
        if x!='-':
            hc+=1
            for k,p in POS.items():
                if p==hc: res[k]=(y.upper() if y!='-' else '-', tc if y!='-' else None)
    # local identity over aligned columns
    pairs=[(x,y) for x,y in zip(h,t) if x!='-' and y!='-']
    ident=sum(1 for x,y in pairs if x.upper()==y.upper())/max(1,len(pairs))
    return acc,res,round(ident,3)
accs=list(M.keys())
with ThreadPoolExecutor(2) as ex: results=list(ex.map(run,accs))
out={acc:dict(res=res,ident=ident) for acc,res,ident in results}
json.dump(out,open('pairmap.json','w'))
print('mapped',len(out))
