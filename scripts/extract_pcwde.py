#!/usr/bin/env python3
# Extract PCWDE-family proteins from dbCAN(HMM) across all species -> per-family FASTA + count matrix
import glob, csv, re, os
from collections import defaultdict
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
DBCAN=PROJ+'/annotation/dbcan'; CLEAN=PROJ+'/clean'
OUTD=PROJ+'/h1_pcwde/families'; os.makedirs(OUTD,exist_ok=True)
ML=PROJ+'/species_lists/manifest_dl.tsv'
PCWDE={'GH5','GH6','GH7','GH10','GH11','GH12','GH28','GH43','GH51','GH53','GH74','GH93',
       'AA9','PL1','PL3','PL4','PL9','CE1','CE5','CE8','CE12','CE15','CE16'}
meta={r['portal']:r for r in csv.DictReader(open(ML),delimiter='\t')}

def fam(h): m=re.match(r'([A-Za-z]+\d+)', h); return m.group(1) if m else None

# protein -> families per species
prot_fam=defaultdict(lambda: defaultdict(set))   # portal -> protein -> {fam}
for hmm in glob.glob(DBCAN+'/*/dbCAN_hmm_results.tsv'):
    portal=hmm.split('/dbcan/')[1].split('/')[0]
    with open(hmm) as f:
        next(f,None)
        for line in f:
            c=line.rstrip('\n').split('\t')
            if len(c)<3: continue
            fm=fam(c[0])
            if fm in PCWDE: prot_fam[portal][c[2]].add(fm)

# load sequences per species (only proteins we need)
def load_seqs(portal, wanted):
    seqs={}; cur=None; buf=[]
    p=CLEAN+'/'+portal+'.faa'
    if not os.path.exists(p): return seqs
    with open(p) as f:
        for line in f:
            if line.startswith('>'):
                if cur in wanted: seqs[cur]=''.join(buf)
                cur=line[1:].strip().split()[0]; buf=[]
            else: buf.append(line.strip())
        if cur in wanted: seqs[cur]=''.join(buf)
    return seqs

fam_records=defaultdict(list)   # family -> [(header, seq)]
count=defaultdict(lambda: defaultdict(int))  # family -> lifestyle -> n
for portal, pmap in prot_fam.items():
    seqs=load_seqs(portal, set(pmap.keys()))
    life=meta.get(portal,{}).get('lifestyle','?')
    for prot, fams in pmap.items():
        s=seqs.get(prot)
        if not s: continue
        for fm in fams:
            fam_records[fm].append(('%s'%prot, s))
            count[fm][life]+=1

for fm, recs in fam_records.items():
    with open(OUTD+'/%s.faa'%fm,'w') as o:
        for h,s in recs: o.write('>%s\n%s\n'%(h,s))

# report
print('PCWDE families extracted:', len(fam_records))
print('%-6s %6s  ECM  SAP  WD  ORM  ERM'%('fam','total'))
for fm in sorted(fam_records, key=lambda f:-len(fam_records[f])):
    c=count[fm]
    print('%-6s %6d  %4d %4d %3d %4d %4d'%(fm,len(fam_records[fm]),c.get('ECM',0),c.get('SAP',0),c.get('WD',0),c.get('ORM',0),c.get('ERM',0)))
