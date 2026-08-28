#!/usr/bin/env python3
# Residue-level catalytic instrument: align family to a characterized reference, read catalytic residues.
# Usage: catalytic_check.py <FAMILY> <ref.fasta> <pos1,pos2,..> <exp1,exp2,..>
import sys, os, subprocess, csv, re
from collections import defaultdict
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
FAM=sys.argv[1]; REFFA=sys.argv[2]
POSITIONS=[int(x) for x in sys.argv[3].split(',')]
EXPECTED=sys.argv[4].split(',')
FAMILIES=PROJ+'/h1_pcwde/families'; MSAD=PROJ+'/h1_pcwde/msa'; os.makedirs(MSAD,exist_ok=True)
ML=PROJ+'/species_lists/manifest_dl.tsv'
meta={r['portal']:r['lifestyle'] for r in csv.DictReader(open(ML),delimiter='\t')}
portals=sorted(meta, key=len, reverse=True)
def portal_of(h):
    for p in portals:
        if h.startswith(p+'_'): return p
    return None

def read_fa(path):
    d={}; cur=None; buf=[]
    for line in open(path):
        if line.startswith('>'):
            if cur: d[cur]=''.join(buf)
            cur=line[1:].strip().split()[0]; buf=[]
        else: buf.append(line.strip())
    if cur: d[cur]=''.join(buf)
    return d

refseq=list(read_fa(REFFA).values())[0]
# build input: REF (POS control) + synthetic N2 (each catalytic mutated) + family
inp=MSAD+'/%s_input.faa'%FAM
with open(inp,'w') as o:
    o.write('>REF_POS\n%s\n'%refseq)
    for i,pos in enumerate(POSITIONS):
        mut=list(refseq); mut[pos-1]='A'
        o.write('>N2_single_%s%d\n%s\n'%(EXPECTED[i],pos,''.join(mut)))
    allmut=list(refseq)
    for pos in POSITIONS: allmut[pos-1]='A'
    o.write('>N2_all\n%s\n'%''.join(allmut))
    for h,s in read_fa(FAMILIES+'/%s.faa'%FAM).items():
        o.write('>%s\n%s\n'%(h,s))
aln=MSAD+'/%s.aln'%FAM
with open(aln,'w') as ao:
    subprocess.run(['/home1/minseo1101/miniforge3/bin/mamba','run','-n','phylogeny','mafft','--auto','--quiet','--thread','16',inp],stdout=ao,stderr=subprocess.DEVNULL)
seqs=read_fa(aln)
refaln=seqs['REF_POS']
# map ref positions -> alignment columns
col={}; rp=0
for c,ch in enumerate(refaln):
    if ch!='-':
        rp+=1
        if rp in POSITIONS: col[rp]=c
cols=[col[p] for p in POSITIONS]
def classify(s):
    res=[s[c] for c in cols]
    if any(r=='-' for r in res): return 'fragment', res
    if all(res[i]==EXPECTED[i] for i in range(len(POSITIONS))): return 'intact', res
    return 'disrupted', res
# validation
print('=== INSTRUMENT VALIDATION (%s: pos %s exp %s) ==='%(FAM,POSITIONS,EXPECTED))
for ctrl in ['REF_POS','N2_single_%s%d'%(EXPECTED[0],POSITIONS[0]),'N2_all']:
    cl,res=classify(seqs[ctrl]); print('  %-24s -> %-9s residues=%s'%(ctrl,cl,res))
# classify family
by=defaultdict(lambda: defaultdict(int)); rows=[]
for h,s in seqs.items():
    if h.startswith('REF_') or h.startswith('N2_'): continue
    cl,res=classify(s); life=meta.get(portal_of(h) or '','?')
    by[life][cl]+=1; rows.append((h,life,cl,''.join(res)))
with open(PROJ+'/h1_pcwde/%s_calls.tsv'%FAM,'w',newline='') as o:
    w=csv.writer(o,delimiter='\t'); w.writerow(['protein','lifestyle','call','catalytic_residues'])
    for r in rows: w.writerow(r)
print('\n=== %s catalytic status by lifestyle ==='%FAM)
print('%-6s %6s %7s %9s %8s  %s'%('life','intact','disrupt','fragment','total','%intact(of non-frag)'))
for life in sorted(by,key=lambda l:-(by[l]['intact']+by[l]['disrupted']+by[l]['fragment'])):
    d=by[life]; intact=d['intact']; dis=d['disrupted']; frag=d['fragment']; tot=intact+dis+frag
    nf=intact+dis; pct=100*intact/nf if nf else 0
    print('%-6s %6d %7d %9d %8d  %.0f%%'%(life,intact,dis,frag,tot,pct))
