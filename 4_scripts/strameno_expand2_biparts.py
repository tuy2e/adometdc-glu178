import sys, json
from Bio import Phylo
def biparts(tree_path):
    t=Phylo.read(tree_path,'newick')
    allt=frozenset(x.name for x in t.get_terminals())
    out={}
    for cl in t.find_clades():
        if cl.is_terminal() or cl==t.root: continue
        s=frozenset(x.name for x in cl.get_terminals())
        key=min(s, allt-s, key=lambda z:(len(z), sorted(z)))
        out[key]=cl.confidence if cl.confidence is not None else (cl.name if cl.name else None)
    return allt,out
def has(allt,bp,S):
    S=frozenset(S); C=allt-S
    for k,v in bp.items():
        if k==S or k==C: return True,v
    return False,None
if __name__=='__main__':
    h=json.load(open('/root/audit/strameno/expand/hyp_sets.json'))
    allt,bp=biparts(sys.argv[1])
    sets={'ST':h['ST'],'DP':h['DP'],'KE':h['KE'],'CR':h['CR']}
    combos={'L=ST+DP+KE+CR':['ST','DP','KE','CR'],'ST+DP':['ST','DP'],'KE+CR':['KE','CR'],'ST':['ST'],'DP':['DP'],'KE':['KE'],'CR':['CR'],'DP+KE+CR':['DP','KE','CR'],'ST+KE+CR':['ST','KE','CR']}
    for n,c in combos.items():
        S=set().union(*[sets[x] for x in c]); S&=allt
        print(n, len(S), has(allt,bp,S))
