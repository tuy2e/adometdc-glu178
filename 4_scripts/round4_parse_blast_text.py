"""Parse NCBI web BLAST 'Download All -> Text' reports into per-hit protein sequences.
For each subject: HSP translations (Sbjct lines, gaps removed). Best HSP plus any non-overlapping HSPs (query coords), joined in query order.
Output: JSON list of records."""
import re, sys, json
def parse(path, search):
    txt=open(path,encoding='utf-8',errors='replace').read()
    meta={}
    for k in ['RID','Job Title','Program','Database']:
        m=re.search(r'^'+re.escape(k)+r':\s*(.*)$',txt,re.M); meta[k]=m.group(1).strip() if m else None
    if 'Alignments:' not in txt: return meta,[]
    body=txt.split('\nAlignments:',1)[1]
    blocks=re.split(r'\n>(?=\S)',('\n'+body))[1:]
    recs=[]
    for b in blocks:
        lines=b.split('\n')
        title=lines[0].strip()
        i=1
        while i<len(lines) and not lines[i].startswith('Sequence ID:'):
            title+=' '+lines[i].strip(); i+=1
        m=re.search(r'Sequence ID:\s*(\S+)\s+Length:\s*(\d+)',b)
        if not m: continue
        acc,slen=m.group(1),int(m.group(2))
        hsps=[]
        for part in re.split(r'\nRange \d+:',b)[1:]:
            sc=re.search(r'Score:\s*([\d.]+) bits',part); ev=re.search(r'Expect:\s*([^,\s]+)',part)
            idm=re.search(r'Identities:\s*(\d+)/(\d+)',part)
            fr=re.search(r'Frame:\s*([+-]?\d+)',part)
            q=[(int(x.group(1)),int(x.group(3))) for x in re.finditer(r'^Query\s+(\d+)\s+(\S+)\s+(\d+)',part,re.M)]
            s=[x.group(2) for x in re.finditer(r'^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)',part,re.M)]
            if not q or not s: continue
            seq=''.join(s).replace('-','').upper()
            hsps.append(dict(bits=float(sc.group(1)) if sc else 0, evalue=ev.group(1) if ev else None,
                             ident=int(idm.group(1))/int(idm.group(2)) if idm else None, alnlen=int(idm.group(2)) if idm else None,
                             qstart=q[0][0], qend=q[-1][1], frame=fr.group(1) if fr else None, seq=seq))
        if not hsps: continue
        hsps.sort(key=lambda h:-h['bits'])
        keep=[hsps[0]]
        for h in hsps[1:]:
            if all(h['qend']<k['qstart']+5 or h['qstart']>k['qend']-5 for k in keep): keep.append(h)
        keep.sort(key=lambda h:h['qstart'])
        sp=None
        mm=re.search(r'\[([^\[\]]+)\]\s*$',title.split(' >')[0])
        if mm: sp=mm.group(1)
        else:
            t=re.sub(r'^(TSA|TPA_asm|TPA):\s*','',title); t=re.sub(r'^uncultured\s+','uncultured_',t)
            w=t.replace(',',' ').split()
            sp=' '.join(w[:2]) if len(w)>1 and w[1][0].islower() else w[0]
        recs.append(dict(search=search, acc=acc, title=title[:200], species=sp, slen=slen, n_hsp=len(hsps),
                         best_bits=hsps[0]['bits'], best_evalue=hsps[0]['evalue'], best_ident=hsps[0]['ident'],
                         qcov=(min(k['qstart'] for k in keep), max(k['qend'] for k in keep)),
                         seq=''.join(k['seq'] for k in keep)))
    return meta,recs
if __name__=='__main__':
    D='/mnt/user-data/uploads/adometdc-glu178/3_data/blast_results/'
    S={'A':'F5G9RPR9016','A2':'F5J4PJTH014','B':'F5J6K5RN014','C':'F5HFU3NK014','D':'F5HK30WR014','D2':'F5J81YRJ014','D3':'F6EP33G2016','E':'F5HNB4Z6014','F':'F5HSK9YB016','G':'F5EJKSE8014','H':'F5FJXSNR014'}
    allr=[]
    for k,rid in S.items():
        meta,recs=parse(D+rid+'-Alignment.txt',k)
        print(k,rid,meta['Job Title'],'|',meta['Program'],'|',len(recs),'hits', 'median len',sorted(len(r['seq']) for r in recs)[len(recs)//2] if recs else 0)
        allr+=recs
    json.dump(allr,open('/home/claude/r4/blast_hits.json','w'),indent=0)
    print('total',len(allr))
