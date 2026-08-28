#!/usr/bin/env python3
# Subsample a (trimmed) codon alignment to N_test ECM + N_ref non-ECM tips for a tractable large-tree RELAX.
# Deterministic (seed=42). Preserves species breadth via random draw within each group.
# Usage: subsample_codon.py <FAM> <N_test> <N_ref>
import sys, random
PROJ='/home1/minseo1101/projects/aicomp-h1h3'; C=PROJ+'/h1_pcwde/codon'
FAM=sys.argv[1]; NT=int(sys.argv[2]); NR=int(sys.argv[3])
random.seed(42)
def read_fa(p):
    d={};cur=None;buf=[]
    for line in open(p):
        if line.startswith('>'):
            if cur:d[cur]=''.join(buf)
            cur=line[1:].strip();buf=[]
        else:buf.append(line.strip())
    if cur:d[cur]=''.join(buf)
    return d
aln=read_fa('%s/%s.trim.codon.fasta'%(C,FAM))
ecm=set(l.strip() for l in open('%s/%s.ecm_tips.txt'%(C,FAM)) if l.strip())
test=[h for h in aln if h in ecm]
ref=[h for h in aln if h not in ecm]
test.sort(); ref.sort()
keep_t = test if len(test)<=NT else random.sample(test,NT)
keep_r = ref if len(ref)<=NR else random.sample(ref,NR)
keep=set(keep_t)|set(keep_r)
with open('%s/%s.sub.codon.fasta'%(C,FAM),'w') as o:
    for h in aln:
        if h in keep: o.write('>%s\n%s\n'%(h,aln[h]))
with open('%s/%s.sub.ecm_tips.txt'%(C,FAM),'w') as o:
    for h in keep_t: o.write(h+'\n')
print('%s subsample: Test %d/%d, Ref %d/%d, total=%d, sites=%d'%(FAM,len(keep_t),len(test),len(keep_r),len(ref),len(keep),len(aln[next(iter(keep))])//3))
