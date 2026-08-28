#!/usr/bin/env python3
# Build the PGLS analysis table: per-species lifestyle, raw PCWDE count, and functional (catalytically intact) count.
import csv, re, os
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
# lifestyle + raw pcwde count + total cazyme
rows={}
for r in csv.DictReader(open(PROJ+'/annotation/pcwde_counts.tsv'),delimiter='\t'):
    rows[r['portal']]={'portal':r['portal'],'lifestyle':r['lifestyle'],
                       'pcwde_count':int(r['pcwde_count']),'total_cazyme':int(r['total_cazyme'])}
# portal resolver from manifest (longest-prefix)
portals=sorted({r['portal'] for r in csv.DictReader(open(PROJ+'/species_lists/manifest_dl.tsv'),delimiter='\t')},key=len,reverse=True)
def pof(h):
    for p in portals:
        if h.startswith(p+'_'): return p
    return None
# functional (intact) counts in the 4 catalytically-assessable families (AA9 excluded: SignalP-dependent, unreliable)
FAMS=['GH5','GH6','GH7','GH28']
func={}; assessed={}
for fam in FAMS:
    for r in csv.DictReader(open(PROJ+'/h1_pcwde/%s_calls.tsv'%fam),delimiter='\t'):
        p=pof(r['protein'])
        if not p: continue
        assessed[p]=assessed.get(p,0)+1
        if r['call']=='intact': func[p]=func.get(p,0)+1
# tree tips
tree=open(PROJ+'/orthofinder/run/Results_Jul14/Species_Tree/SpeciesTree_rooted.txt').read()
tips=set(re.findall(r'[\(,]([A-Za-z0-9_]+):',tree))
# assemble, restricted to species present in the tree
out=[]
for p,d in rows.items():
    if p not in tips: continue
    d['func_count']=func.get(p,0)
    d['assessed_count']=assessed.get(p,0)
    out.append(d)
cols=['portal','lifestyle','pcwde_count','total_cazyme','func_count','assessed_count']
with open(PROJ+'/h3_pgls/pgls_data.tsv','w') as o:
    o.write('\t'.join(cols)+'\n')
    for d in sorted(out,key=lambda x:x['portal']):
        o.write('\t'.join(str(d[c]) for c in cols)+'\n')
print('species in tree ∩ counts:', len(out))
# quick lifestyle summary
from collections import defaultdict
agg=defaultdict(lambda:[0,0,0])
for d in out:
    a=agg[d['lifestyle']]; a[0]+=1; a[1]+=d['pcwde_count']; a[2]+=d['func_count']
print('%-8s %4s %10s %10s'%('lifestyle','n','pcwde_mean','func_mean'))
for k in sorted(agg,key=lambda k:-agg[k][1]/agg[k][0]):
    n,s,f=agg[k]; print('%-8s %4d %10.1f %10.1f'%(k,n,s/n,f/n))
print('missing tree tips (in counts, not tree):', sorted(set(rows)-tips)[:10])
