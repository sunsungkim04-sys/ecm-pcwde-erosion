#!/usr/bin/env python3
# Request restoration of PURGED CDS files from JGI tape archive (request_archived_files endpoint).
import sys, json, subprocess, os
sys.path.insert(0,'/home1/minseo1101/projects/aicomp-h1h3/scripts')
import resolve_portals as R
PROJ='/home1/minseo1101/projects/aicomp-h1h3'
TOKEN=open('/home1/minseo1101/.jgi_token').read().strip()
passed=[l.strip() for l in open(PROJ+'/qc/qc_pass.txt') if l.strip()]
ids={}; have=0; nocds=0
for portal in passed:
    out=PROJ+'/cds/%s.cds.fna.gz'%portal
    if os.path.exists(out) and subprocess.run(['gzip','-t',out],capture_output=True).returncode==0:
        have+=1; continue
    d=R.search(portal)
    if not d: continue
    for o in d.get('organisms',[]):
        if R.portal_of(o)!=portal: continue
        th=(o.get('top_hit') or {}).get('_id')
        cds=[f for f in o.get('files',[]) if 'cds' in (f.get('file_type') or [])
             and 'genecatalog' in f.get('file_name','').lower() and f.get('file_name','').endswith('.fasta.gz')]
        if not cds:
            cds=[f for f in o.get('files',[]) if 'cds' in (f.get('file_type') or []) and f.get('file_name','').endswith('.fasta.gz')]
        if cds:
            ids[o.get('id')]={'file_ids':[cds[0]['_id']],'top_hit':th,'mycocosm_portal_id':portal}
        else: nocds+=1
        break
body={'ids':ids,'send_mail':False,'api_version':'2'}
bf=PROJ+'/logs/restore_body.json'; open(bf,'w').write(json.dumps(body))
print('requesting restore for %d files (already have %d, no_cds %d)'%(len(ids),have,nocds))
r=subprocess.run(['wget','--header=Authorization: '+TOKEN,'--header=Content-Type: application/json',
                  '--post-file='+bf,'--content-on-error','-q','-O','-',
                  'https://files.jgi.doe.gov/request_archived_files/'],capture_output=True,timeout=300)
resp=r.stdout.decode()[:600]
print('RESPONSE:', resp)
open(PROJ+'/logs/restore_response.json','w').write(r.stdout.decode())
