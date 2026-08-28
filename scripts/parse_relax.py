#!/usr/bin/env python3
# Parse HyPhy RELAX JSON -> k, p-value, LRT, and Test/Reference omega distributions.
# Usage: parse_relax.py <fam>.relax.json
import sys, json
p=sys.argv[1]
d=json.load(open(p))
fam=p.split('/')[-1].split('.')[0]

tr=d.get('test results',{})
k=tr.get('relaxation or intensification parameter')
pval=tr.get('p-value')
lrt=tr.get('LRT')

def omega(model):
    # returns list of (omega, weight) from a fitted model's Test/Reference rate distribution
    fits=d.get('fits',{}).get(model,{})
    rd=fits.get('Rate Distributions',{})
    out={}
    for branchset in rd:  # e.g. 'Test','Reference'
        cats=rd[branchset]
        pts=[]
        for kk in sorted(cats, key=lambda x:int(x) if str(x).isdigit() else x):
            c=cats[kk]
            pts.append((round(c.get('omega',float('nan')),4), round(c.get('proportion',float('nan')),4)))
        out[branchset]=pts
    return out

alt=omega('RELAX alternative')
# branch counts
bl=d.get('tested',{}) or d.get('branch attributes',{})
print(f"=== {fam} ===")
print(f"k (relaxation/intensification) = {k}")
print(f"p-value = {pval}   LRT = {lrt}")
if k is not None:
    if pval is not None and pval<0.05:
        verdict = "RELAXED selection on Test (k<1)" if k<1 else "INTENSIFIED selection on Test (k>1)"
    else:
        verdict = "no significant difference (function maintained on Test branches)"
    print(f"interpretation = {verdict}")
print("omega distributions (RELAX alternative):")
for bs,pts in alt.items():
    print(f"  {bs}: " + ", ".join(f"w={w}(p={pr})" for w,pr in pts))
# input tip / branch summary
inp=d.get('input',{})
print(f"input: sequences={inp.get('number of sequences')} sites={inp.get('number of sites')}")
