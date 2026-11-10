"""Posterior states at named nodes for chosen human positions. usage: node_states.py prefix colmap.json node=label ..."""
import sys, json
pre=sys.argv[1]; cm={int(k):v for k,v in json.load(open(sys.argv[2])).items()}
nodes=dict(x.split('=') for x in sys.argv[3:])
POS=[13,15,109,111,174,176,178,256]
cols={cm[p]:p for p in POS}
hdr=None; out={}
for l in open(pre+'.state'):
    if l.startswith('#') or not l.strip(): continue
    f=l.rstrip('\n').split('\t')
    if f[0]=='Node': hdr=f; continue
    if f[0] in nodes and int(f[1]) in cols:
        pr=sorted(((float(f[i]),hdr[i][2:]) for i in range(3,len(hdr))),reverse=True)
        out.setdefault(f[0],{})[cols[int(f[1])]]=' '.join('%s %.2f'%(a,p) for p,a in pr[:3] if p>=0.05)
for n,lab in nodes.items():
    print('%-8s %-34s'%(n,lab)+' | '.join('%d: %s'%(p,out.get(n,{}).get(p,'?')) for p in POS))
