import json
# lineages for AF-derived entries (transcribed from UniProt taxonomy query, 5 Nov)
L={'2903':'Haptista, Haptophyta','3075':'Viridiplantae, Chlorophyta','5693':'Discoba, Euglenozoa, Kinetoplastea, Trypanosomatida',
'7757':'Opisthokonta, Metazoa, Chordata','13221':'Haptista, Haptophyta','31354':'Rhodophyta','46634':'Opisthokonta, Fungi',
'53468':'Opisthokonta, Metazoa, Platyhelminthes, Cestoda','59798':'Discoba, Euglenozoa, Kinetoplastea, Trypanosomatida',
'63605':'Discoba, Heterolobosea','72004':'Opisthokonta, Metazoa, Chordata','73025':'Discoba, Euglenozoa, Euglenida',
'89044':'Sar, Stramenopiles, Ochrophyta, Chrysophyceae','94130':'Opisthokonta, Fungi','124430':'Sar, Stramenopiles, Ochrophyta, Pinguiophyceae',
'156173':'Haptista, Haptophyta','156174':'Haptista, Haptophyta','216820':'Sar, Stramenopiles, Ochrophyta, Bacillariophyta, Fragilariophyceae',
'312471':'Discoba, Euglenozoa, Kinetoplastea, Neobodonida','420281':'Sar, Stramenopiles, Ochrophyta, Bacillariophyta, Melosirophyceae',
'464990':'Cryptophyceae','582737':'Viridiplantae, Chlorophyta','641309':'Sar, Rhizaria, Cercozoa','652834':'Eukaryota incertae sedis, Palpitomonas',
'765719':'Viridiplantae, Chlorophyta','861768':'Opisthokonta, Metazoa, Chordata','867726':'Sar, Stramenopiles, Ochrophyta, Phaeophyceae',
'1169539':'Sar, Alveolata, Colpodellida','1411642':'Viridiplantae, Chlorophyta','1561963':'Sar, Rhizaria, Cercozoa','1622150':'Opisthokonta, Fungi',
'1714386':'Sar, Stramenopiles, Ochrophyta, Bolidophyceae','2483409':'Sar, Stramenopiles, Oomycota','2768841':'Opisthokonta, Fungi',
'3111310':'Sar, Stramenopiles, Ochrophyta, Dictyochophyceae'}
master={}
def add(e,org,lin,seq,src,frag='',pfam=''):
    if e in master: return
    master[e]=dict(acc=e,org=org,lineage=lin,seq=seq,len=len(seq),src=src,fragment=frag,pfam=pfam)
for f,src in [('uniprot_strameno.json','UniProt'),('sar_other.json','UniProt'),('outgroups.json','UniProt')]:
    for r in json.load(open(f)):
        add(r['Entry'],r['Organism'],r['Taxonomic lineage'],r['Sequence'],src,r['Fragment'],r['Pfam'])
for r in json.load(open('af_derived.json')):
    add(r['Entry'],r['Organism'],L[r['taxid']],r['Sequence'],'AFmodel')
def classify(lin):
    l=lin
    if 'Stramenopiles' in l:
        if 'Ochrophyta' in l: return ('Stramenopiles','Ochrophyta')
        if 'Oomycota' in l: return ('Stramenopiles','Oomycota')
        if 'Bigyra' in l: return ('Stramenopiles','Bigyra')
        return ('Stramenopiles','other')
    if 'Alveolata' in l:
        for k in ['Ciliophora','Apicomplexa','Dinophyceae','Perkinsozoa','Colpodellida','Chromerida']:
            if k in l: return ('Alveolata',k)
        return ('Alveolata','other')
    if 'Rhizaria' in l: return ('Rhizaria','Rhizaria')
    if 'Haptophyta' in l: return ('Haptophyta','Haptophyta')
    if 'Cryptophyceae' in l: return ('Cryptophyta','Cryptophyceae')
    if 'Palpitomonas' in l: return ('Cryptista-like','Palpitomonas')
    if 'Rhodophyta' in l: return ('Archaeplastida','Rhodophyta')
    if 'Viridiplantae' in l:
        if 'Embryophyta' in l: return ('Archaeplastida','land plants')
        if 'Chlorophyta' in l: return ('Archaeplastida','Chlorophyta')
        return ('Archaeplastida','other green')
    if 'Euglenozoa' in l:
        if 'Kinetoplastea' in l: return ('Discoba','Kinetoplastea')
        return ('Discoba','Euglenida')
    if 'Heterolobosea' in l: return ('Discoba','Heterolobosea')
    if 'Metazoa' in l: return ('Opisthokonta','Metazoa')
    if 'Fungi' in l: return ('Opisthokonta','Fungi')
    if 'Opisthokonta' in l: return ('Opisthokonta','other opisthokonts')
    if 'Amoebozoa' in l or 'Evosea' in l or 'Discosea' in l: return ('Amoebozoa','Amoebozoa')
    if 'Parabasalia' in l or 'Metamonada' in l: return ('Metamonada','Parabasalia')
    return ('?','?')
from collections import Counter
for m in master.values(): m['group'],m['sub']=classify(m['lineage'])
print('total sequences',len(master))
print(Counter(m['group'] for m in master.values()))
print('unclassified:',[(m['acc'],m['org'],m['lineage'][:80]) for m in master.values() if m['group']=='?'])
json.dump(master,open('master.json','w'))
with open('master.fa','w') as f:
    for m in master.values(): f.write('>%s\n%s\n'%(m['acc'],m['seq']))
