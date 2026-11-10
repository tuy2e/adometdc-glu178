# Add Egr_SpdSyn to the original SpdSyn AU-test constraint trees.
# Euglena belongs with kinetoplastids in AdoMetDC group (b), so it goes inside any constrained group that contains kinetoplastids.
import re, json
W='/home/claude/ws/audit/spdsyn/'; O='/home/claude/eug/spd/'
EGR='Egr_SpdSyn_GDJR01022181'
lab=json.load(open(W+'labels2.json'))
KE={k for k,v in lab.items() if v['sub']=='Kinetoplastea'}
for h in ['h1_crke','h2_crrho','h3_stdino','h4_lineage']:
    s=open(W+h+'.nwk').read().strip()
    assert s.startswith('((') and s.endswith(');'), h
    inner_end=s.index(')')
    group=s[2:inner_end].split(','); rest=s[inner_end+1:-2].strip(',').split(',')
    assert '(' not in s[2:inner_end] and s.count('(')==2, 'unexpected nesting in '+h
    if KE & set(group): group.append(EGR); where='inside constrained group'
    else: rest.append(EGR); where='unconstrained'
    open(O+h+'.nwk','w').write('(('+','.join(group)+'),'+','.join(rest)+');\n')
    print(h, len(group),'constrained +',len(rest),'free; Egr',where)
