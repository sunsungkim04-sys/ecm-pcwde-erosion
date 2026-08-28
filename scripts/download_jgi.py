#!/usr/bin/env python3
# Download proteome + assembly per species from JGI Data Portal; verify md5 + gzip; idempotent.
import json, subprocess, csv, os, zipfile, hashlib, sys, time
PROJ='/home1/minseo1101/projects/aicomp-h1h3/'
BASE=PROJ+'species_lists/'; PROT=PROJ+'proteomes/'; GEN=PROJ+'genomes/'; TMP=PROJ+'tmp_zips/'; LOG=PROJ+'logs/'
TOKEN=open('/home1/minseo1101/.jgi_token').read().strip()
for d in (PROT,GEN,TMP,LOG): os.makedirs(d,exist_ok=True)

def md5(path):
    h=hashlib.md5()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''): h.update(chunk)
    return h.hexdigest()

def gzip_ok(path):
    return subprocess.run(['gzip','-t',path],capture_output=True).returncode==0

def post_download(org_id, fids, zippath, tries=3):
    body=json.dumps({"ids":{org_id:fids},"api_version":"2"})
    for i in range(tries):
        r=subprocess.run(['wget','--header=Authorization: '+TOKEN,'--header=Content-Type: application/json',
                          '--post-data='+body,'--content-on-error','-q','-O',zippath,
                          'https://files-download.jgi.doe.gov/download_files/'],capture_output=True,timeout=600)
        if os.path.exists(zippath) and zipfile.is_zipfile(zippath): return True
        time.sleep(2)
    return False

def extract_member(zf, wanted_name, out):
    for n in zf.namelist():
        if n.split('/')[-1]==wanted_name:
            with zf.open(n) as src, open(out,'wb') as dst: dst.write(src.read())
            return True
    return False

def process(row):
    portal=row['portal']; org=row['org_id']
    pf,af=row['proteome_fid'],row['assembly_fid']
    pmd5,amd5=row.get('proteome_md5',''),row.get('assembly_md5','')
    prot_out=PROT+portal+'.faa.gz'; gen_out=GEN+portal+'.fna.gz'
    # idempotent skip
    def good(path,exp):
        return os.path.exists(path) and gzip_ok(path) and (not exp or md5(path)==exp)
    if good(prot_out,pmd5) and good(gen_out,amd5): return 'SKIP'
    zp=TMP+portal+'.zip'
    if not post_download(org,[pf,af],zp): return 'DL_FAIL'
    try:
        zf=zipfile.ZipFile(zp)
        ok1=extract_member(zf,row['proteome_name'],prot_out)
        ok2=extract_member(zf,row['assembly_name'],gen_out)
    except Exception as e:
        return 'ZIP_FAIL:%s'%e
    finally:
        try: os.remove(zp)
        except: pass
    if not (ok1 and ok2): return 'MISSING_MEMBER'
    if not (gzip_ok(prot_out) and gzip_ok(gen_out)): return 'GZIP_BAD'
    if pmd5 and md5(prot_out)!=pmd5: return 'MD5_PROT_MISMATCH'
    if amd5 and md5(gen_out)!=amd5: return 'MD5_ASM_MISMATCH'
    return 'OK'

if __name__=='__main__':
    rows=[r for r in csv.DictReader(open(BASE+"manifest_dl.tsv"),delimiter='\t') if r.get('status')=='OK']
    log=open(LOG+'download.log','a')
    counts={}; 
    for i,r in enumerate(rows,1):
        st=process(r); counts[st.split(':')[0]]=counts.get(st.split(':')[0],0)+1
        log.write('%s\t%s\t%s\n'%(r['portal'],r['name'],st)); log.flush()
        if st not in ('OK','SKIP'): print('  !!',r['portal'],r['name'],st,flush=True)
        if i%20==0: print('...%d/%d %s'%(i,len(rows),counts),flush=True)
    print('DONE',counts)
