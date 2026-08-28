#!/usr/bin/env python3
# Trim gap-heavy codon columns from a codon alignment to make large RELAX runs tractable.
# Keeps every sequence; drops codon columns whose gap fraction exceeds --maxgap (default 0.80).
# Usage: trim_codon.py <in.codon.fasta> <out.codon.fasta> [maxgap]
import sys
inp, outp = sys.argv[1], sys.argv[2]
maxgap = float(sys.argv[3]) if len(sys.argv) > 3 else 0.80
def read_fa(p):
    d={};cur=None;buf=[]
    for line in open(p):
        if line.startswith('>'):
            if cur:d[cur]=''.join(buf)
            cur=line[1:].strip();buf=[]
        else:buf.append(line.strip())
    if cur:d[cur]=''.join(buf)
    return d
d=read_fa(inp)
names=list(d)
L=len(d[names[0]])
assert L%3==0, "alignment length not multiple of 3"
ncod=L//3
n=len(names)
keep=[]
for c in range(ncod):
    col=[d[h][c*3:c*3+3] for h in names]
    gaps=sum(1 for x in col if x=='---')
    if gaps/n <= maxgap: keep.append(c)
with open(outp,'w') as o:
    for h in names:
        s=d[h]
        o.write('>%s\n%s\n'%(h,''.join(s[c*3:c*3+3] for c in keep)))
print(f"{inp}: seqs={n} codons {ncod}->{len(keep)} (maxgap={maxgap}) new_len={len(keep)*3}")
