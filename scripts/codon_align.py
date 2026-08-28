#!/usr/bin/env python3
# Back-translate a family protein MSA to a codon alignment, restricted to catalytically-intact proteins.
# Usage: codon_align.py <FAMILY>
#
# CDS<->protein mapping fix: JGI assigns proteinId to the proteome but transcriptId to the CDS FASTA;
# the two numeric IDs differ per release (only the shared model name, 4th '|' field, is stable).
# We therefore bridge proteinId -> model name (from proteomes/) -> nt (from cds/), and fall back to
# translated-sequence matching when a model name is absent/unshared.
import sys, os, glob, gzip, csv, re
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
FAM=sys.argv[1]
MSAD=PROJ+'/h1_pcwde/msa'; OUTD=PROJ+'/h1_pcwde/codon'; os.makedirs(OUTD,exist_ok=True)
CDSDIR=PROJ+'/cds'; PROTD=PROJ+'/proteomes'
ML=PROJ+'/species_lists/manifest_dl.tsv'
meta={r['portal']:r['lifestyle'] for r in csv.DictReader(open(ML),delimiter='\t')}
portals=sorted(meta,key=len,reverse=True)
def portal_of(h):
    for p in portals:
        if h.startswith(p+'_'): return p
    return None

def read_fa(path, gz=False):
    op=gzip.open(path,'rt') if gz else open(path)
    d={}; cur=None; buf=[]
    for line in op:
        if line.startswith('>'):
            if cur: d[cur]=''.join(buf)
            cur=line[1:].strip().split()[0]; buf=[]
        else: buf.append(line.strip())
    if cur: d[cur]=''.join(buf)
    return d

def parse_jgi(h):
    # 'jgi|Portal|ID|modelname...' -> (ID, modelname); modelname None if absent
    parts=h.split('|')
    if len(parts)>=4 and parts[0]=='jgi': return parts[2], parts[3]
    if len(parts)>=3 and parts[0]=='jgi': return parts[2], None
    return h, None

_GC={
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def translate(nt):
    nt=nt.upper()
    return ''.join(_GC.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3))

# intact proteins for this family
intact=set()
for r in csv.DictReader(open(PROJ+'/h1_pcwde/%s_calls.tsv'%FAM),delimiter='\t'):
    if r['call']=='intact': intact.add(r['protein'])

# alignment (protein), keep only intact real proteins
aln=read_fa(MSAD+'/%s.aln'%FAM)
prots={h:s for h,s in aln.items() if h in intact}
print('%s: intact aligned proteins=%d'%(FAM,len(prots)))

# ungapped family protein sequences (for translated-sequence fallback)
protseq={h:aap.replace('-','').upper() for h,aap in prots.items()}

# build CDS map keyed by {portal}_{proteinId}, bridging proteinId<->transcriptId via JGI model name
need_portals=set(filter(None,(portal_of(h) for h in prots)))
cdsmap={}
via_model=0; via_pep=0
for portal in need_portals:
    fp=CDSDIR+'/%s.cds.fna.gz'%portal
    if not os.path.exists(fp): continue
    cds=read_fa(fp, gz=True)
    model2nt={}; pep2nt={}
    for h,s in cds.items():
        _, model = parse_jgi(h)
        if model: model2nt.setdefault(model, s)
        pep=translate(s).rstrip('*')
        pep2nt.setdefault(pep, s)
    # proteinId -> model name from proteome (full JGI headers)
    prot2model={}
    ppath=PROTD+'/%s.faa.gz'%portal
    if os.path.exists(ppath):
        for h in read_fa(ppath, gz=True):
            pid, model = parse_jgi(h)
            if model: prot2model[pid]=model
    plen=len(portal)+1
    for cleanid, seq in protseq.items():
        if portal_of(cleanid)!=portal: continue
        pid=cleanid[plen:]
        nt=None
        model=prot2model.get(pid)
        if model and model in model2nt: nt=model2nt[model]
        if nt is not None:
            cdsmap[cleanid]=nt; via_model+=1; continue
        nt=pep2nt.get(seq)                       # fallback: exact translated-peptide match
        if nt is not None:
            cdsmap[cleanid]=nt; via_pep+=1
print('%s: cds mapped via model=%d via translation=%d'%(FAM,via_model,via_pep))

def backtranslate(aap, nt):
    # strip trailing stop
    if len(nt)>=3 and nt[-3:].upper() in ('TAA','TAG','TGA'): nt=nt[:-3]
    codons=[nt[i:i+3] for i in range(0,len(nt),3)]
    out=[]; ci=0
    for aa in aap:
        if aa=='-': out.append('---')
        else:
            if ci>=len(codons): return None
            out.append(codons[ci]); ci+=1
    if ci!=len(codons): return None   # length mismatch = model/CDS inconsistent
    return ''.join(out)

recs=[]; miss=0; mism=0
for h,aap in prots.items():
    nt=cdsmap.get(h)
    if not nt: miss+=1; continue
    cod=backtranslate(aap, nt)
    if cod is None: mism+=1; continue
    recs.append((h,cod))
with open(OUTD+'/%s.codon.fasta'%FAM,'w') as o:
    for h,s in recs: o.write('>%s\n%s\n'%(h,s))
# labels for RELAX (ECM tips = Test)
with open(OUTD+'/%s.ecm_tips.txt'%FAM,'w') as o:
    for h,_ in recs:
        if meta.get(portal_of(h))=='ECM': o.write(h+'\n')
print('%s: codon-aligned=%d  (no_cds=%d, len_mismatch=%d) ECM=%d'%(FAM,len(recs),miss,mism,sum(1 for h,_ in recs if meta.get(portal_of(h))=='ECM')))
