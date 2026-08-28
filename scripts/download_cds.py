#!/usr/bin/env python3
# Resolve + download CDS (nucleotide) for QC-passed species -> cds/<portal>.cds.fna.gz (for RELAX codon alignment)
import json, subprocess, csv, os, zipfile, hashlib, sys, re, time
sys.path.insert(0,'/home1/minseo1101/projects/aicomp-h1h3/scripts')
import resolve_portals as R
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
CDS=PROJ+'/cds'; TMP=PROJ+'/tmp_zips'; LOG=PROJ+'/logs'
os.makedirs(CDS,exist_ok=True); os.makedirs(TMP,exist_ok=True)
TOKEN=open('/home1/minseo1101/.jgi_token').read().strip()
passed=[l.strip() for l in open(PROJ+'/qc/qc_pass.txt') if l.strip()]
meta={r['portal']:r for r in csv.DictReader(open(PROJ+'/species_lists/manifest_dl.tsv'),delimiter='\t')}

def pick_cds(files):
    cands=[]
    for f in files:
        nl=f.get('file_name','').lower(); ft=f.get('file_type') or []
        if 'cds' not in ft or not nl.endswith('.fasta.gz'): continue
        # prefer GeneCatalog CDS / FilteredModels CDS, avoid transcripts/ESTs
        if 'transcript' in nl or 'est' in nl: continue
        score = (2 if 'cds' in nl else 0) + (1 if 'genecatalog' in nl or 'filteredmodels' in nl else 0)
        cands.append((score, f.get('file_size',0) or 0, f))
    if not cands: return None
    cands.sort(key=lambda x:(x[0],x[1]), reverse=True); return cands[0][2]

def gzip_ok(p): return subprocess.run(['gzip','-t',p],capture_output=True).returncode==0
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for c in iter(lambda:f.read(1<<20),b''): h.update(c)
    return h.hexdigest()

def resolve_cds(portal):
    d=R.search(portal)
    if not d: return None
    for org in d.get('organisms',[]):
        if R.portal_of(org)!=portal: continue
        cds=pick_cds(org.get('files',[]))
        if cds: return org.get('id'), cds
    # fallback: any org with cds file matching portal in filename
    for org in d.get('organisms',[]):
        cds=pick_cds(org.get('files',[]))
        if cds and cds.get('file_name','').lower().startswith(portal.lower()):
            return org.get('id'), cds
    return None

log=open(LOG+'/download_cds.log','a'); counts={}
for i,portal in enumerate(passed,1):
    out=CDS+'/'+portal+'.cds.fna.gz'
    if os.path.exists(out) and gzip_ok(out): counts['SKIP']=counts.get('SKIP',0)+1; continue
    r=resolve_cds(portal)
    if not r: counts['NO_CDS']=counts.get('NO_CDS',0)+1; log.write('%s\tNO_CDS\n'%portal); continue
    org,cds=r; body=json.dumps({"ids":{org:[cds['_id']]},"api_version":"2"})
    zp=TMP+'/'+portal+'_cds.zip'; ok=False
    for _ in range(3):
        subprocess.run(['wget','--header=Authorization: '+TOKEN,'--header=Content-Type: application/json','--post-data='+body,'--content-on-error','-q','-O',zp,'https://files-download.jgi.doe.gov/download_files/'],capture_output=True,timeout=600)
        if os.path.exists(zp) and zipfile.is_zipfile(zp): ok=True; break
        time.sleep(2)
    if not ok: counts['DL_FAIL']=counts.get('DL_FAIL',0)+1; log.write('%s\tDL_FAIL\n'%portal); continue
    zf=zipfile.ZipFile(zp); wrote=False
    for n in zf.namelist():
        if n.split('/')[-1]==cds['file_name']:
            open(out,'wb').write(zf.open(n).read()); wrote=True; break
    os.remove(zp)
    if wrote and gzip_ok(out) and (not cds.get('md5sum') or md5(out)==cds['md5sum']):
        counts['OK']=counts.get('OK',0)+1
    else:
        counts['VERIFY_FAIL']=counts.get('VERIFY_FAIL',0)+1; log.write('%s\tVERIFY_FAIL\n'%portal)
    if i%20==0: print('...%d/%d %s'%(i,len(passed),counts),flush=True)
print('CDS DONE',counts)
