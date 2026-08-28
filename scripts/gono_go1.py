#!/usr/bin/env python3
# Go/No-Go 1: PCWDE raw counts per species -> ECM vs free-living decomposer (Miyauchi reduction check)
import glob, csv, re, os
from collections import defaultdict
try:
    from scipy.stats import mannwhitneyu
    HAVE_SCIPY=True
except Exception:
    HAVE_SCIPY=False
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
DBCAN=PROJ+'/annotation/dbcan'
ML=PROJ+'/species_lists/manifest_dl.tsv'
PCWDE={'GH5','GH6','GH7','GH10','GH11','GH12','GH28','GH43','GH51','GH53','GH74','GH93',
       'AA9','PL1','PL3','PL4','PL9','CE1','CE5','CE8','CE12','CE15','CE16'}
meta={r['portal']:r for r in csv.DictReader(open(ML),delimiter='\t')}

def fam(hmmname):
    m=re.match(r'([A-Za-z]+\d+)', hmmname); return m.group(1) if m else None

counts={}         # portal -> {family: set(proteins)}
for hmm in glob.glob(DBCAN+'/*/dbCAN_hmm_results.tsv'):
    portal=hmm.split('/dbcan/')[1].split('/')[0]
    fp=defaultdict(set)
    with open(hmm) as f:
        next(f,None)
        for line in f:
            c=line.rstrip('\n').split('\t')
            if len(c)<3: continue
            fm=fam(c[0])
            if fm: fp[fm].add(c[2])
    counts[portal]=fp

rows=[]
for p,fp in counts.items():
    life=meta.get(p,{}).get('lifestyle','?')
    pcwde=sum(len(fp.get(f,())) for f in PCWDE)
    total_caz=len(set().union(*fp.values())) if fp else 0
    rows.append((p,life,pcwde,total_caz))
# write per-species
with open(PROJ+'/annotation/pcwde_counts.tsv','w',newline='') as o:
    w=csv.writer(o,delimiter='\t'); w.writerow(['portal','lifestyle','pcwde_count','total_cazyme'])
    for r in sorted(rows,key=lambda x:x[1]): w.writerow(r)

from statistics import mean, median
grp=defaultdict(list)
for p,life,pcwde,_ in rows: grp[life].append(pcwde)
print('=== PCWDE raw count by lifestyle (n, mean, median) ===')
for life in sorted(grp, key=lambda l:-mean(grp[l])):
    v=grp[life]; print('  %-6s n=%-3d mean=%5.1f median=%4.1f'%(life,len(v),mean(v),median(v)))

ecm=grp.get('ECM',[]); free=grp.get('SAP',[])+grp.get('WD',[])
print('\n=== Go/No-Go 1 ===')
print('  species assessed (dbCAN done): %d'%len(rows))
print('  ECM: %d  | free-living decomposer (SAP+WD): %d'%(len(ecm),len(free)))
if ecm and free:
    print('  ECM mean PCWDE=%.1f vs free-living=%.1f'%(mean(ecm),mean(free)))
    if HAVE_SCIPY:
        u,pval=mannwhitneyu(ecm,free,alternative='less')
        print('  Mann-Whitney (ECM < free-living): U=%.0f  p=%.2e  -> Miyauchi reduction %s'%(u,pval,'REPRODUCED' if pval<0.05 else 'NOT sig'))
