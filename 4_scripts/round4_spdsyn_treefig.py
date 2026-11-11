"""Collapsed cladogram of the spermidine synthase tree, rooted on prokaryotes.
Uniform clades collapse to one labelled row; colours mark the groups that matter for the AdoMetDC comparison."""
import json, sys, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
tf, labf, out, title = sys.argv[1:5]
L=json.load(open(labf)); t=Phylo.read(tf,'newick')
G=lambda a:L[a]['grp']; SB=lambda a:L[a]['sub']
def key(a):
    g,s=G(a),SB(a)
    if g=='DISC': return {'Kinetoplastea':'kinetoplastids','Euglenida':'euglenids','Diplonemea':'diplonemids','Heterolobosea':'heteroloboseans'}.get(s,s)
    if g=='CRY': return 'Goniomonas (colourless cryptophyte)' if L[a]['org'].startswith(('Goniomonas','uncultured Goniomonas')) else 'cryptophytes'
    if g=='ALV': return {'Dinophyceae':'dinoflagellates','Perkinsozoa':'Perkinsus','Apicomplexa':'apicomplexans','Colpodellida':'chromerids','Chromerida':'chromerids'}.get(s,'alveolates')
    if g=='ST': return 'oomycetes' if s=='Oomycota' else ('diatoms' if s=='Bacillariophyta' else ('Bigyra' if s=='Bigyra' else 'other stramenopiles'))
    return {'KAT':'katablepharid (Roombia)','HAP':'haptophytes','RHI':'Rhizaria','RHO':'red algae','VIR':'green algae + plants','OPI':'animals + fungi','AMO':'amoebozoans','MET':'metamonads','PROK':'prokaryotes'}.get(g,g)
COL={'kinetoplastids':'#c2410c','cryptophytes':'#1d4ed8','euglenids':'#b45309','diplonemids':'#7c3aed','Goniomonas (colourless cryptophyte)':'#60a5fa',
     'dinoflagellates':'#047857','Perkinsus':'#047857','diatoms':'#0f766e','oomycetes':'#0f766e','Bigyra':'#0f766e','other stramenopiles':'#0f766e'}
P=[n for n in t.get_terminals() if G(n.name)=='PROK']; ca=t.common_ancestor(P)
t.root_with_outgroup(ca if {x.name for x in ca.get_terminals()}=={n.name for n in P} else P[0])
t.ladderize()
rows=[]  # (y, x_end, label, color, n)
segs=[]; y=[0]
def sup(c):
    v=c.name if c.name else (str(c.confidence) if c.confidence is not None else '')
    try: return float(v.split('/')[-1])
    except: return None
def draw(c,x0):
    tips=[t_.name for t_ in c.get_terminals()]
    ks=collections.Counter(key(a) for a in tips)
    x1=x0+min(c.branch_length or 0,1.2)
    if len(ks)==1:
        k=list(ks)[0]; yy=y[0]; y[0]+=1
        rows.append((yy,x1,f"{k} ({len(tips)})" if len(tips)>1 else f"{k}: {L[tips[0]]['org'][:28]}",COL.get(k,'#6b7280'),len(tips)))
        segs.append(((x0,yy),(x1,yy),COL.get(k,'#9ca3af'))); return yy
    ys=[draw(ch,x1) for ch in c.clades]
    yy=(min(ys)+max(ys))/2
    segs.append(((x1,min(ys)),(x1,max(ys)),'#9ca3af')); segs.append(((x0,yy),(x1,yy),'#9ca3af'))
    s=sup(c)
    if s is not None and s>=90 and len(tips)>=4: rows.append((yy,x1,'__sup__%d'%s,'#111827',0))
    return yy
draw(t.root,0)
n=y[0]; fig,ax=plt.subplots(figsize=(10,max(6,n*0.22)))
for (a,b,c) in segs: ax.plot([a[0],b[0]],[a[1],b[1]],color=c,lw=1.2,solid_capstyle='butt')
xmax=max(r[1] for r in rows)
for yy,x,lab,col,nn in rows:
    if lab.startswith('__sup__'): ax.text(x-0.02,yy-0.15,lab[7:],fontsize=6,color='#374151',ha='right',va='bottom'); continue
    ax.text(x+0.03,yy,lab,fontsize=7.5,va='center',color=col,fontweight='bold' if col!='#6b7280' else 'normal')
ax.set_ylim(n,-1); ax.set_xlim(0,xmax+1.6); ax.axis('off')
ax.set_title(title,fontsize=10,loc='left')
ax.text(0,n+0.5,'Rooted on prokaryotes. Numbers = UFBoot support (≥90 shown). Uniform clades are collapsed; counts in brackets. Branch lengths capped at 1.2 for display.',fontsize=7,color='#4b5563')
plt.tight_layout(); plt.savefig(out,dpi=200)
print('rows',n)
