#!/usr/bin/env python3
# Consolidate all RELAX results into one table.
# Usage: summarize_relax.py   (scans h1_pcwde/codon/*.relax.json)
import json, glob, os
C='/home1/minseo1101/projects/aicomp-h1h3/h1_pcwde/codon'
rows=[]
for jp in sorted(glob.glob(C+'/*.relax.json')):
    if os.path.getsize(jp)==0: continue
    d=json.load(open(jp))
    name=os.path.basename(jp).replace('.relax.json','')
    tr=d.get('test results',{})
    k=tr.get('relaxation or intensification parameter')
    pv=tr.get('p-value'); lrt=tr.get('LRT')
    inp=d.get('input',{})
    nseq=inp.get('number of sequences'); nsit=inp.get('number of sites')
    # count Test branches from branch attributes / tested
    ntest=None
    tested=d.get('tested',{}).get('0',{})
    if tested:
        ntest=sum(1 for v in tested.values() if str(v).lower()=='test')
    if k is None:
        rows.append((name,nseq,ntest,nsit,'NA','NA','NA','no result')); continue
    if pv is not None and pv<0.05:
        verd='RELAXED (k<1)' if k<1 else 'INTENSIFIED (k>1)'
    else:
        verd='maintained (ns)'
    rows.append((name,nseq,ntest,nsit,round(k,3),('%.2g'%pv if pv is not None else 'NA'),(round(lrt,2) if lrt is not None else 'NA'),verd))
print('%-14s %6s %6s %6s %7s %10s %8s  %s'%('analysis','nseq','nTest','sites','k','p','LRT','verdict'))
for r in rows:
    print('%-14s %6s %6s %6s %7s %10s %8s  %s'%r)
