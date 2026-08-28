#!/usr/bin/env python3
# Annotate ECM tips with {Test} in the trimmed-alignment Newick tree for HyPhy RELAX.
import sys, re
PROJ='/home1/minseo1101/projects/aicomp-h1h3'; FAM=sys.argv[1]
C=PROJ+'/h1_pcwde/codon'
tree=open('%s/%s.trim.tree.nwk'%(C,FAM)).read()
ecm=set(l.strip() for l in open('%s/%s.ecm_tips.txt'%(C,FAM)) if l.strip())
def repl(m):
    name=m.group(1)
    return name+('{Test}' if name in ecm else '')+':'
tree2=re.sub(r'(?<=[(,])([A-Za-z0-9_.\-]+):', repl, tree)
open('%s/%s.trim.labeled.nwk'%(C,FAM),'w').write(tree2)
print('%s: labeled %d ECM tips as {Test}'%(FAM, tree2.count('{Test}')))
