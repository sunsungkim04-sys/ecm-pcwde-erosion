#!/usr/bin/env python3
# Re-collate BUSCO + apply filter: PASS = BUSCO complete >=85% (proteome-based; N50 informational).
import csv, os, glob, re
QC='/home1/minseo1101/projects/aicomp-h1h3/qc'
ML='/home1/minseo1101/projects/aicomp-h1h3/species_lists/manifest_dl.tsv'
# re-collate busco
busco={}
for f in glob.glob(QC+'/busco/*/short_summary*.txt'):
    portal=f.split('/busco/')[1].split('/')[0]
    m=re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]',open(f).read())
    if m: busco[portal]=float(m.group(1))
astat={}
for r in csv.DictReader(open(QC+'/assembly_stats.tsv'),delimiter='\t'):
    p=os.path.basename(r['file']).replace('.fna.gz',''); astat[p]=(int(r['N50']),int(r['num_seqs']))
meta={r['portal']:r for r in csv.DictReader(open(ML),delimiter='\t')}
rows=[]; passed=[]
for p in sorted(set(list(busco)+list(astat))):
    c=busco.get(p); n50,ct=astat.get(p,(0,10**9))
    ok = c is not None and c>=85.0
    frag = not((n50>=100000) or (ct<5000))
    rows.append((p, meta.get(p,{}).get('lifestyle','?'), c if c is not None else -1, n50, ct, 'frag' if frag else '', 'PASS' if ok else 'FAIL'))
    if ok: passed.append(p)
with open(QC+'/qc_final.tsv','w',newline='') as o:
    w=csv.writer(o,delimiter='\t'); w.writerow(['portal','lifestyle','busco_complete','N50','contigs','flag','result'])
    for r in rows: w.writerow(r)
open(QC+'/qc_pass.txt','w').write('\n'.join(passed))
from collections import Counter
print('assessed:',len(rows),'| PASS(BUSCO>=85):',len(passed),'| FAIL:',len(rows)-len(passed))
print('PASS by lifestyle:', dict(Counter(meta.get(p,{}).get('lifestyle','?') for p in passed)))
print('  ECM pass:', sum(1 for p in passed if meta.get(p,{}).get('lifestyle')=='ECM'))
print('FAILS:', [(r[0],r[2]) for r in rows if r[6]=='FAIL'])
print('fragmented-but-passed (proteome OK, assembly frag):', sum(1 for r in rows if r[5]=='frag' and r[6]=='PASS'))
