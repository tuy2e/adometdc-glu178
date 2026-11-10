"""Tree figure: ML tree of active AdoMetDCs, branches coloured by the reconstructed residue at
human position 178; uniform clades collapsed to wedges. usage: python3 treefig.py <treefile> <statefile> <site> <out.png> <title>"""
import sys, json, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from Bio import Phylo

treef, statef, site, outpng, title = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
M = json.load(open('master_annot.json'))
alnf = sys.argv[6] if len(sys.argv) > 6 else 'tree_trim.aln'
aln = {}; k = None
for l in open(alnf):
    l = l.strip()
    if l.startswith('>'): k = l[1:].split()[0]; aln[k] = ''
    else: aln[k] += l
anc = {}; hdr = None
if statef != '-':
    for l in open(statef):
        if l.startswith('#') or not l.strip(): continue
        f = l.rstrip('\n').split('\t')
        if f[0] == 'Node': hdr = f; continue
        if int(f[1]) == site: anc[f[0]] = (f[2], max(float(x) for x in f[3:]))
# tip residue = multiple-alignment column at the human-178 position
tipres = {a: aln[a][site - 1].upper() for a in aln}  # same column the ASR used

tree = Phylo.read(treef, 'newick')
for cl in tree.get_nonterminals():
    p = (cl.name or '').split('/')
    cl.node = p[0] if p and p[0].startswith('Node') else None
    cl.ufb = float(p[2]) if len(p) > 2 and p[2] not in ('', 'nan') else None
if getattr(tree.root, 'node', None) is None and 'Node1' in anc: tree.root.node = 'Node1'
# root: between the opisthokont+amoebozoan block and the rest is not available (not monophyletic),
# so root on human P17707's clade side using midpoint for display only
tree.root_at_midpoint()
tree.ladderize(reverse=True)

COL = {'E': '#b9b8b2', 'N': '#2a78d6', 'H': '#eb6834', 'D': '#1baf7a'}
OTHER = '#52514e'
def col(s): return '#d6d5cf' if s == '?' else COL.get(s, OTHER)
def st(cl):
    if cl.is_terminal(): return tipres[cl.name]
    return anc.get(getattr(cl, 'node', None), ('?', 0))[0]

SHORT = {'Ochrophyta': 'Ochrophyta', 'Oomycota': 'Oomycota', 'Bigyra': 'Bigyra', 'Metazoa': 'animals',
         'Fungi': 'fungi', 'land plants': 'land plants', 'Chlorophyta': 'green algae', 'Rhodophyta': 'red algae',
         'Haptophyta': 'haptophytes', 'Dinophyceae': 'dinoflagellates', 'Kinetoplastea': 'kinetoplastids',
         'Euglenida': 'euglenid', 'Cryptophyceae': 'cryptophyte', 'Heterolobosea': 'heteroloboseans',
         'Perkinsozoa': 'Perkinsus', 'Apicomplexa': 'apicomplexans', 'Colpodellida': 'chromerids',
         'Rhizaria': 'rhizarians', 'Amoebozoa': 'amoebozoans', 'Parabasalia': 'parabasalid',
         'other opisthokonts': 'other opisthokonts', 'other green': 'other green'}
STRAM = {'Ochrophyta', 'Oomycota', 'Bigyra'}

def uniform(cl):
    t = cl.get_terminals()
    subs = {M[x.name]['sub'] for x in t}; res = {tipres[x.name] for x in t}
    return len(subs) == 1 and len(res) == 1

# layout: collapse uniform clades with >=3 tips into wedges
rows = []   # (clade, collapsed?)
def layout(cl):
    if cl.is_terminal() or (len(cl.get_terminals()) >= 3 and uniform(cl)):
        rows.append(cl); return
    for c in cl.clades: layout(c)
layout(tree.root)
depth = tree.depths()
y = {}; ROWH = {}
cur = 0
for cl in rows:
    h = 1.0 if cl.is_terminal() else 1.0 + 0.18 * min(len(cl.get_terminals()), 30)
    y[cl] = cur + h / 2; ROWH[cl] = h; cur += h
def sety(cl):
    if cl in y: return y[cl]
    ys = [sety(c) for c in cl.clades]; y[cl] = (min(ys) + max(ys)) / 2; return y[cl]
sety(tree.root)
H = cur
fig, ax = plt.subplots(figsize=(11, max(8, H * 0.16)))
xmax = max(depth.values())
changes = []
def draw(cl):
    x0 = depth[cl]
    if cl in rows and not cl.is_terminal():
        t = cl.get_terminals(); s = tipres[t[0].name]
        far = max(depth[x] for x in t)
        h = ROWH[cl] * 0.42
        ax.add_patch(Polygon([[x0, y[cl]], [far, y[cl] - h], [far, y[cl] + h]], closed=True,
                             facecolor=col(s), edgecolor='white', lw=1.5, alpha=0.9))
        sub = M[t[0].name]['sub']
        lab = f"{SHORT.get(sub, sub)} ({len(t)})"
        ax.text(far + xmax * 0.01, y[cl], lab, va='center', fontsize=8.5,
                fontweight='bold' if sub in STRAM else 'normal', color='#0b0b0b')
        ax.text(xmax * 1.33, y[cl], s, va='center', ha='center', fontsize=9, fontweight='bold', color='#0b0b0b')
        return
    if cl.is_terminal():
        s = tipres[cl.name]; sub = M[cl.name]['sub']
        org = M[cl.name]['org'].split('(')[0].strip()
        org = ' '.join(org.split()[:2])
        ax.plot([x0], [y[cl]], 'o', ms=4.5, color=col(s), mec='white', mew=0.8, zorder=3)
        ax.text(x0 + xmax * 0.01, y[cl], f"{org}  [{SHORT.get(sub, sub)}]", va='center', fontsize=7.2,
                style='italic', fontweight='bold' if sub in STRAM else 'normal', color='#0b0b0b')
        ax.text(xmax * 1.33, y[cl], s, va='center', ha='center', fontsize=8, color='#0b0b0b')
        return
    ys = [y[c] for c in cl.clades]
    ax.plot([x0, x0], [min(ys), max(ys)], color=col(st(cl)), lw=1.6, solid_capstyle='round')
    for c in cl.clades:
        ax.plot([x0, depth[c]], [y[c], y[c]], color=col(st(c)), lw=1.6, solid_capstyle='round')
        ps, cs = st(cl), st(c)
        if ps != cs and '?' not in (ps, cs) and '-' not in (ps, cs) and ps == 'E':
            changes.append((x0 + (depth[c] - x0) * 0.5, y[c], ps, cs, getattr(c, 'ufb', None)))
        draw(c)
    if getattr(cl, 'ufb', None) is not None and cl.ufb >= 95 and len(cl.get_terminals()) >= 6:
        ax.text(x0 - xmax * 0.004, y[cl] + 0.35, f"{cl.ufb:.0f}", fontsize=6, ha='right', va='bottom', color='#52514e')
draw(tree.root)
for (cx, cy, a, b, u) in changes:
    big = u is not None
    ax.plot([cx], [cy], marker='D', ms=9 if big else 6, color='#0b0b0b', mec='white', mew=1, zorder=5)
    if big:
        ax.text(cx, cy - 0.6, f"Glu→Asn on this branch\n(UFBoot {u:.0f}; 63 without pocket sites)", fontsize=8.5, ha='center', va='bottom', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85))
    else:
        ax.text(cx, cy + 0.25, f"{a}→{b}", fontsize=7, ha='center', va='top', color='#0b0b0b',
                bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.85))
ax.set_ylim(H + 1, -1)
ax.set_xlim(-xmax * 0.02, xmax * 1.36)
ax.text(xmax * 1.33, -0.6, 'res.\n178', ha='center', va='bottom', fontsize=8, color='#52514e')
ax.axis('off')
ax.text(xmax * 0.27 + xmax, -0.6, 'res. 178', ha='center', fontsize=8, color='#52514e') if False else None
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], color=COL['E'], lw=3, label='Glu (E) — the putrescine sensor'),
           Line2D([0], [0], color=COL['N'], lw=3, label='Asn (N)'),
           Line2D([0], [0], color=COL['H'], lw=3, label='His (H)'),
           Line2D([0], [0], color=COL['D'], lw=3, label='Asp (D) — still acidic'),
           Line2D([0], [0], color=OTHER, lw=3, label='other (S, Q, R, A …)'),
           Line2D([0], [0], marker='D', color='#0b0b0b', lw=0, ms=7, label='reconstructed loss of Glu178')]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.0, 0.985), fontsize=8.5, frameon=False, title='Residue at human position 178\n(branches = ancestral reconstruction)', title_fontsize=8.5)
ax.set_title(title, fontsize=11, loc='left')
ax.plot([0, 0.2], [H + 0.3, H + 0.3], color='#0b0b0b', lw=1)
ax.text(0.1, H + 0.6, '0.2 subst./site', ha='center', va='top', fontsize=7)
plt.tight_layout()
plt.savefig(outpng, dpi=200)
print('rows', len(rows), 'height', H, 'changes', [(a, b, u) for (_, _, a, b, u) in changes])
