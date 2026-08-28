#!/usr/bin/env python3
# Sanitize JGI proteome headers -> >{portal}_{proteinId}; write to clean/<portal>.faa
import gzip, os, glob, re
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
SRC=PROJ+'/proteomes'; DST=PROJ+'/clean'; os.makedirs(DST,exist_ok=True)
n=0
for f in sorted(glob.glob(SRC+'/*.faa.gz')):
    portal=os.path.basename(f).replace('.faa.gz','')
    out=DST+'/'+portal+'.faa'
    if os.path.exists(out) and os.path.getsize(out)>0: continue
    seen=set(); i=0
    with gzip.open(f,'rt') as inp, open(out,'w') as o:
        for line in inp:
            if line.startswith('>'):
                parts=line[1:].strip().split('|')
                if len(parts)>=3 and parts[0]=='jgi':
                    pid=re.sub(r'[^A-Za-z0-9]','_',parts[2].split()[0])
                else:
                    pid=re.sub(r'[^A-Za-z0-9]','_',line[1:].split()[0]) if line[1:].strip() else str(i)
                tag='%s_%s'%(portal,pid)
                while tag in seen: i+=1; tag='%s_%s_%d'%(portal,pid,i)
                seen.add(tag); o.write('>'+tag+'\n')
            else:
                o.write(line)
    n+=1
print('cleaned proteomes written:', n, '(existing skipped)')
print('total in clean/:', len(glob.glob(DST+'/*.faa')))
