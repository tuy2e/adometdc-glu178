import json,sys
def parse(path):
    d=json.load(open(path)); txt=''.join(x['text'] for x in d)
    body=txt[txt.index('---\n')+4:]
    trunc = '\nTab Context:' not in body
    if not trunc: body=body[:body.index('\nTab Context:')]
    lines=[l for l in body.split('\n') if l.strip()]
    hdr=lines[0].split('\t'); rows=[]; bad=[]
    for l in lines[1:]:
        p=l.split('\t')
        if len(p)!=len(hdr): bad.append(l[:60]); continue
        r=dict(zip(hdr,p))
        if len(r['Sequence'])!=int(r['Length']): bad.append(r['Entry']+' len-mismatch'); continue
        rows.append(r)
    return rows,bad,len(txt)
if __name__=='__main__':
    rows,bad,n=parse(sys.argv[1])
    print('chars',n,'| rows ok',len(rows),'| bad',bad)
    if len(sys.argv)>2: json.dump(rows,open(sys.argv[2],'w'))
