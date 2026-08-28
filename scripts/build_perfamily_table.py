#!/usr/bin/env python3
# Per-species, per-family intact PCWDE counts (for the H1<->H3 bridge PGLS).
import csv, re
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
portals=sorted({r['portal'] for r in csv.DictReader(open(PROJ+'/species_lists/manifest_dl.tsv'),delimiter='\t')},key=len,reverse=True)
def pof(h):
    for p in portals:
        if h.startswith(p+'_'): return p
    return None
life={r['portal']:r['lifestyle'] for r in csv.DictReader(open(PROJ+'/annotation/pcwde_counts.tsv'),delimiter='\t')}
FAMS=['GH7','GH6','GH28','GH5']
# per portal per family intact count
cnt={p:{f:0 for f in FAMS} for p in life}
for fam in FAMS:
    for r in csv.DictReader(open(PROJ+'/h1_pcwde/%s_calls.tsv'%fam),delimiter='\t'):
        if r['call']!='intact': continue
        p=pof(r['protein'])
        if p in cnt: cnt[p][fam]+=1
# tree tips
tree=open(PROJ+'/orthofinder/run/Results_Jul14/Species_Tree/SpeciesTree_rooted.txt').read()
tips=set(re.findall(r'[\(,]([A-Za-z0-9_]+):',tree))
cols=['portal','lifestyle']+FAMS
with open(PROJ+'/h3_pgls/perfamily_intact.tsv','w') as o:
    o.write('\t'.join(cols)+'\n')
    for p in sorted(cnt):
        if p not in tips: continue
        o.write('\t'.join([p,life[p]]+[str(cnt[p][f]) for f in FAMS])+'\n')
# summary: ECM vs SAP+WD mean per family + % species with >=1 intact
import statistics as st
grp={'ECM':[],'DEC':[]}
for p in cnt:
    if p not in tips: continue
    g='ECM' if life[p]=='ECM' else ('DEC' if life[p] in ('SAP','WD') else None)
    if g: grp[g].append(cnt[p])
print('%-5s %8s %8s %8s %8s'%('fam','ECM_mean','DEC_mean','ECM_%pres','DEC_%pres'))
for f in FAMS:
    em=st.mean(c[f] for c in grp['ECM']); dm=st.mean(c[f] for c in grp['DEC'])
    ep=100*sum(1 for c in grp['ECM'] if c[f]>0)/len(grp['ECM'])
    dp=100*sum(1 for c in grp['DEC'] if c[f]>0)/len(grp['DEC'])
    print('%-5s %8.2f %8.2f %7.0f%% %7.0f%%'%(f,em,dm,ep,dp))
