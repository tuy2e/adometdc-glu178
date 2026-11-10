"""Ancestral-state analysis of the human-Glu178 column on the IQ-TREE ML tree.
Inputs: tree178.treefile (node labels written by -asr), tree178.state, tree_trim.aln,
master_annot.json, site_cols.json. Output: asr178.json + printed summary."""
import json, re, collections, sys
PRE = sys.argv[1] if len(sys.argv) > 1 else 'tree178'
from Bio import Phylo

SITE = json.load(open('site_cols.json'))['E178']          # 1-based column in tree_trim.aln
M = json.load(open('master_annot.json'))

# --- tip residues straight from the alignment column (independent of pairmap) ---
aln = {}; k = None
for l in open('tree_trim.aln'):
    l = l.strip()
    if l.startswith('>'): k = l[1:].split()[0]; aln[k] = ''
    else: aln[k] += l
tipres = {a: s[SITE - 1].upper() for a, s in aln.items()}
assert tipres['P17707'] == 'E', 'column check failed: human is not E at site %d' % SITE

# --- ancestral states at that site ---
anc = {}
hdr = None
for l in open(PRE + '.state'):
    if l.startswith('#') or not l.strip(): continue
    f = l.rstrip('\n').split('\t')
    if f[0] == 'Node': hdr = f; continue
    if int(f[1]) != SITE: continue
    probs = {hdr[i][2:]: float(f[i]) for i in range(3, len(hdr))}
    anc[f[0]] = (f[2], probs[f[2]] if f[2] in probs else max(probs.values()), probs)

tree = Phylo.read(PRE + '.treefile', 'newick')
# internal labels look like NodeN/SH-aLRT/UFBoot (or NodeN only)
for cl in tree.get_nonterminals():
    lab = cl.name or (str(cl.confidence) if cl.confidence is not None else '')
    parts = lab.split('/')
    cl.node = parts[0] if parts and parts[0].startswith('Node') else None
    cl.sh = float(parts[1]) if len(parts) > 1 and parts[1] not in ('', 'nan') else None
    cl.ufb = float(parts[2]) if len(parts) > 2 and parts[2] not in ('', 'nan') else None

if getattr(tree.root, 'node', None) is None and 'Node1' in anc:
    tree.root.node = 'Node1'; tree.root.sh = tree.root.ufb = None
grp = {a: (M[a]['group'], M[a]['sub']) for a in tipres}
def tips(cl): return [t.name for t in cl.get_terminals()]

# root on Amorphea (Opisthokonta + Amoebozoa) if monophyletic, else Opisthokonta
amorph = [a for a in tipres if grp[a][0] in ('Opisthokonta', 'Amoebozoa')]
# first root on an arbitrary non-Amorphea tip so the Amorphea test is meaningful
tree.root_with_outgroup('A0A1E7FKD8')
def is_clade(sel):
    ca = tree.common_ancestor(*sel); return set(tips(ca)) == set(sel), ca
ok, ca = is_clade(amorph)
if ok:
    tree.root_with_outgroup(ca); rooted_on = 'Amorphea (monophyletic in ML tree)'
else:
    opis = [a for a in tipres if grp[a][0] == 'Opisthokonta']
    ok2, ca2 = is_clade(opis)
    tree.root_with_outgroup(ca2 if ok2 else opis[0])
    rooted_on = 'Opisthokonta' + (' (monophyletic)' if ok2 else ' NOT monophyletic; rooted on single tip')
    print('Amorphea not monophyletic; intruders:', sorted(set(tips(ca)) - set(amorph))[:10])

def state(cl):
    if cl.is_terminal(): return tipres[cl.name], 1.0
    s = anc.get(getattr(cl, "node", None))
    return (s[0], s[1]) if s else ('?', 0)

# --- monophyly / support of key groups ---
def clade_of(sel):
    ca = tree.common_ancestor(*sel)
    members = set(tips(ca))
    return ca, members == set(sel), sorted(members - set(sel))
groups = {
 'Stramenopiles': [a for a in tipres if grp[a][0] == 'Stramenopiles'],
 'Ochrophyta': [a for a in tipres if grp[a][1] == 'Ochrophyta'],
 'Oomycota': [a for a in tipres if grp[a][1] == 'Oomycota'],
 'Bigyra': [a for a in tipres if grp[a][1] == 'Bigyra'],
 'Ochrophyta+Bigyra': [a for a in tipres if grp[a][1] in ('Ochrophyta', 'Bigyra')],
 'Kinetoplastea': [a for a in tipres if grp[a][1] == 'Kinetoplastea'],
 'Alveolata': [a for a in tipres if grp[a][0] == 'Alveolata'],
 'Dinophyceae': [a for a in tipres if grp[a][1] == 'Dinophyceae'],
}
res = {'rooted_on': rooted_on, 'site': SITE, 'groups': {}, 'changes': []}
for g, sel in groups.items():
    if len(sel) < 2: continue
    ca, mono, intruders = clade_of(sel)
    st, p = state(ca)
    path = tree.get_path(ca)
    parent = path[-2] if len(path) > 1 else tree.root
    pst, pp = state(parent)
    res['groups'][g] = dict(n=len(sel), monophyletic=mono, intruders=[(a, grp[a][1], tipres[a]) for a in intruders],
                            node=ca.node, sh=ca.sh, ufb=ca.ufb, mrca_state=st, mrca_p=round(p, 3),
                            parent_state=pst, parent_p=round(pp, 3),
                            tip_states=dict(collections.Counter(tipres[a] for a in sel)))

# --- count state changes along every edge (MAP states) ---
for parent in tree.find_clades(order='preorder'):
    if parent.is_terminal(): continue
    ps, pp = state(parent)
    for ch in parent.clades:
        cs, cp = state(ch)
        if ps != cs and '-' not in (ps, cs) and '?' not in (ps, cs):
            t = tips(ch)
            comp = collections.Counter(grp[a][1] for a in t)
            res['changes'].append(dict(frm=ps, frm_p=round(pp, 3), to=cs, to_p=round(cp, 3), n_tips=len(t),
                                       clade=dict(comp), example=t[:3], child_ufb=getattr(ch, 'ufb', None) if not ch.is_terminal() else 'tip'))
acid = set('DE')
res['E_to_nonE'] = [c for c in res['changes'] if c['frm'] == 'E']
res['acid_to_nonacid'] = [c for c in res['changes'] if c['frm'] in acid and c['to'] not in acid]
res['nonE_to_E'] = [c for c in res['changes'] if c['to'] == 'E']
Phylo.write(tree, PRE + '_rooted.nwk', 'newick')
json.dump(res, open(PRE + '_asr.json', 'w'), indent=1)

print('rooted on', rooted_on, '| site', SITE)
for g, d in res['groups'].items():
    print(f"{g:18s} n={d['n']:3d} mono={d['monophyletic']} SH={d['sh']} UFB={d['ufb']} MRCA={d['mrca_state']}({d['mrca_p']}) parent={d['parent_state']}({d['parent_p']}) tips={d['tip_states']} intruders={d['intruders'][:5]}")
print('\nAll MAP state changes on edges:')
for c in res['changes']:
    print(f"  {c['frm']}({c['frm_p']}) -> {c['to']}({c['to_p']})  tips={c['n_tips']:3d} UFB={c['child_ufb']}  {c['clade']}")
print('\nE -> non-E edges:', len(res['E_to_nonE']), '| acidic -> non-acidic edges:', len(res['acid_to_nonacid']), '| -> E edges:', len(res['nonE_to_E']))
