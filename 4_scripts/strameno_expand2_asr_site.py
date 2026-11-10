"""ASR summary at the human-Glu178 column for the expanded trees.
usage: python3 asr_exp.py <prefix> <alignment> [site]"""
import json, sys, collections
from Bio import Phylo
PRE, ALN = sys.argv[1], sys.argv[2]
SITE = int(sys.argv[3]) if len(sys.argv) > 3 else json.load(open('exp_site_cols.json'))['E178']
S = json.load(open('exp_set.json'))
aln = {}; k = None
for l in open(ALN):
    l = l.strip()
    if l.startswith('>'): k = l[1:]; aln[k] = ''
    else: aln[k] += l
tipres = {a: s[SITE - 1].upper() for a, s in aln.items()}
#assert
anc = {}; hdr = None
for l in open(PRE + '.state'):
    if l.startswith('#') or not l.strip(): continue
    f = l.rstrip('\n').split('\t')
    if f[0] == 'Node': hdr = f; continue
    if int(f[1]) != SITE: continue
    probs = {hdr[i][2:]: float(f[i]) for i in range(3, len(hdr))}
    anc[f[0]] = (f[2], probs.get(f[2], 0.0))
tree = Phylo.read(PRE + '.treefile', 'newick')
for cl in tree.get_nonterminals():
    p = (cl.name or '').split('/')
    cl.node = p[0] if p and p[0].startswith('Node') else None
    cl.sh = float(p[1]) if len(p) > 1 and p[1] not in ('', 'nan') else None
    cl.ufb = float(p[2]) if len(p) > 2 and p[2] not in ('', 'nan') else None
if getattr(tree.root, 'node', None) is None and 'Node1' in anc: tree.root.node = 'Node1'
G = lambda a: S[a]['group']; SB = lambda a: S[a]['sub']
tips = lambda c: [t.name for t in c.get_terminals()]
T = list(tipres)
# root
nonamorph = [a for a in T if G(a) not in ('Opisthokonta', 'Amoebozoa', 'Apusozoa')][0]
tree.root_with_outgroup(nonamorph)
def isclade(sel):
    ca = tree.common_ancestor(*sel); return set(tips(ca)) == set(sel), ca
rooted = None
for name, sel in [('Amorphea', [a for a in T if G(a) in ('Opisthokonta', 'Amoebozoa', 'Apusozoa')]),
                  ('Opisthokonta', [a for a in T if G(a) == 'Opisthokonta'])]:
    ok, ca = isclade(sel)
    if ok: tree.root_with_outgroup(ca); rooted = name + ' (monophyletic)'; break
if not rooted: tree.root_with_outgroup('P17707'); rooted = 'human tip (Amorphea/Opisthokonta not monophyletic)'
def st(c):
    if c.is_terminal(): return tipres[c.name], 1.0
    return anc.get(getattr(c, 'node', None), ('?', 0))
groups = collections.OrderedDict([
 ('Stramenopiles', [a for a in T if G(a) == 'Stramenopiles']),
 ('Ochrophyta', [a for a in T if SB(a) == 'Ochrophyta']),
 ('Oomycota', [a for a in T if SB(a) == 'Oomycota']),
 ('Bigyra/other', [a for a in T if G(a) == 'Stramenopiles' and SB(a) not in ('Ochrophyta', 'Oomycota')]),
 ('Cryptophyceae', [a for a in T if SB(a) == 'Cryptophyceae']),
 ('Kinetoplastea', [a for a in T if SB(a) == 'Kinetoplastea']),
 ('Euglenozoa', [a for a in T if SB(a) in ('Kinetoplastea', 'Euglenida', 'Diplonemea')]),
 ('Dinophyceae (standalone)', [a for a in T if SB(a) == 'Dinophyceae']),
 ('Dino+Perkinsus', [a for a in T if SB(a) in ('Dinophyceae', 'Perkinsozoa')]),
 ('Dinophyceae (ODC fusion)', [a for a in T if SB(a) == 'Dinophyceae (ODC fusion)']),
 ('Alveolata', [a for a in T if G(a) == 'Alveolata']),
 ('Haptophyta', [a for a in T if G(a) == 'Haptophyta']),
 ('Rhizaria', [a for a in T if G(a) == 'Rhizaria']),
])
res = {'rooted': rooted, 'site': SITE, 'groups': {}, 'changes': []}
for g, sel in groups.items():
    sel = [a for a in sel if a in tipres]
    if len(sel) < 2: continue
    ca = tree.common_ancestor(*sel); mem = set(tips(ca))
    s, p = st(ca); par = tree.get_path(ca); par = par[-2] if len(par) > 1 else tree.root
    ps, pp = st(par)
    res['groups'][g] = dict(n=len(sel), mono=mem == set(sel), n_clade=len(mem),
                            intruders=dict(collections.Counter(SB(a) for a in mem - set(sel))),
                            sh=getattr(ca, 'sh', None), ufb=getattr(ca, 'ufb', None),
                            mrca=f'{s}({p:.3f})', parent=f'{ps}({pp:.3f})',
                            tips=dict(collections.Counter(tipres[a] for a in sel)))
for par in tree.find_clades(order='preorder'):
    if par.is_terminal(): continue
    ps, pp = st(par)
    for ch in par.clades:
        cs, cp = st(ch)
        if ps != cs and '-' not in (ps, cs) and '?' not in (ps, cs):
            t = tips(ch)
            res['changes'].append(dict(frm=ps, frm_p=round(pp, 3), to=cs, to_p=round(cp, 3), n=len(t),
                                       ufb=getattr(ch, 'ufb', None) if not ch.is_terminal() else 'tip',
                                       clade=dict(collections.Counter(SB(a) for a in t)),
                                       ex=[S[a]['org'][:25] for a in t[:2]]))
json.dump(res, open(PRE + '_asr.json', 'w'), indent=1)
Phylo.write(tree, PRE + '_rooted.nwk', 'newick')
print('rooted on', rooted)
for g, d in res['groups'].items():
    print(f"{g:26s} n={d['n']:3d} mono={str(d['mono']):5s} clade={d['n_clade']:3d} SH={d['sh']} UFB={d['ufb']} MRCA={d['mrca']} parent={d['parent']} tips={d['tips']} intruders={d['intruders']}")
print('\nstate changes (MAP):')
for c in res['changes']:
    print(f"  {c['frm']}({c['frm_p']}) -> {c['to']}({c['to_p']}) n={c['n']:3d} UFB={c['ufb']} {c['clade']} {c['ex'] if c['n'] < 3 else ''}")
E = [c for c in res['changes'] if c['frm'] == 'E']; toE = [c for c in res['changes'] if c['to'] == 'E']
print('\nE->non-E:', len(E), ' ->E:', len(toE), ' acidic->non-acidic:', len([c for c in res['changes'] if c['frm'] in 'DE' and c['to'] not in 'DE']))
