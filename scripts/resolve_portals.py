#!/usr/bin/env python3
# Resolve each species -> JGI org_id + proteome/assembly file_id via Data Portal search.
import json, subprocess, urllib.parse, re, sys, csv, time
BASE='/home1/minseo1101/projects/aicomp-h1h3/species_lists/'
TOKEN=open('/home1/minseo1101/.jgi_token').read().strip()

def search(q, tries=4):
    url='https://files.jgi.doe.gov/search/?q=%s&x=25&api_version=2'%urllib.parse.quote(q)
    for i in range(tries):
        try:
            r=subprocess.run(['wget','--header=Authorization: '+TOKEN,'--content-on-error','-q','-O','-',url],
                             capture_output=True, timeout=90)
            d=json.loads(r.stdout)
            if isinstance(d,dict) and 'organisms' in d: return d
        except Exception:
            pass
    return None

def norm(s):
    s=str(s).lower(); s=re.sub(r'\bv[\d.]+\b',' ',s); s=re.sub(r'[^a-z0-9 ]',' ',s)
    return re.sub(r'\s+',' ',s).strip()

def portal_of(org):
    m=(org.get('top_hit',{}) or {}).get('metadata',{}) or {}
    if m.get('mycocosm_portal_id'): return m['mycocosm_portal_id']
    return None

def pick_proteome(files):
    cands=[]
    for f in files:
        nl=f.get('file_name','').lower(); ft=f.get('file_type') or []
        if 'protein' not in ft or not nl.endswith('.fasta.gz') or 'protein' not in nl: continue
        if any(b in nl for b in ['deflines','secondary','allmodels','_all_','repeat']): continue
        score=0
        if 'genecatalog_proteins' in nl: score+=3
        if 'filteredmodels' in nl and 'proteins' in nl: score+=3
        if nl.endswith('.aa.fasta.gz'): score+=1
        cands.append((score, f.get('file_size',0) or 0, f))
    if not cands: return None
    cands.sort(key=lambda x:(x[0],x[1]), reverse=True); return cands[0][2]

def pick_assembly(files):
    cands=[]
    for f in files:
        nl=f.get('file_name','').lower(); ft=f.get('file_type') or []
        if 'assembly' not in ft or not nl.endswith('.fasta.gz'): continue
        score=0
        if 'assemblyscaffolds' in nl: score+=2
        if 'masked' in nl or 'repeatmask' in nl: score-=1
        cands.append((score, f.get('file_size',0) or 0, f))
    if not cands: return None
    cands.sort(key=lambda x:(x[0],x[1]), reverse=True); return cands[0][2]

def resolve_query(query, want_portal=None, want_name=None):
    d=search(query)
    if not d: return None
    wn=norm(want_name) if want_name else None
    scored=[]
    for org in d.get('organisms',[]):
        files=org.get('files',[])
        prot=pick_proteome(files); asm=pick_assembly(files)
        if not (prot and asm): continue
        pid=portal_of(org); s=0
        if want_portal and pid==want_portal: s+=100
        if wn and norm(org.get('name',''))==wn: s+=50
        elif wn and wn in norm(org.get('name','')): s+=20
        scored.append((s,org,pid,prot,asm))
    if not scored: return None
    scored.sort(key=lambda x:x[0], reverse=True)
    s,org,pid,prot,asm=scored[0]
    return dict(portal=pid, org_id=org.get('id'), org_name=org.get('name'),
                proteome_fid=prot.get('_id'), proteome_name=prot.get('file_name'),
                proteome_md5=prot.get('md5sum'), proteome_size=prot.get('file_size'),
                assembly_fid=asm.get('_id'), assembly_name=asm.get('file_name'),
                assembly_md5=asm.get('md5sum'), assembly_size=asm.get('file_size'), match_score=s)

def resolve_species(name, portal):
    # try portal code first (fast+exact), then full name, then genus+species
    queries=[]
    if portal: queries.append((portal, portal))
    queries.append((name, portal or None))
    toks=norm(name).split()
    if len(toks)>=2: queries.append((' '.join(toks[:2]), portal or None))
    for q,wp in queries:
        res=resolve_query(q, want_portal=wp, want_name=name)
        if res: return res
    return None

if __name__=='__main__':
    rows=list(csv.DictReader(open(BASE+'manifest_pre.tsv'), delimiter='\t'))
    fields=['name','lifestyle','ecology','order','published','source','portal','org_id',
            'proteome_fid','proteome_name','proteome_md5','proteome_size',
            'assembly_fid','assembly_name','assembly_md5','assembly_size','status']
    out=open(BASE+'manifest.tsv','w',newline=''); w=csv.DictWriter(out,fieldnames=fields,delimiter='\t'); w.writeheader()
    unresolved=[]; ok=0
    for i,r in enumerate(rows,1):
        name=r['name']; portal=r.get('portal','').strip()
        res=resolve_species(name, portal or None)
        rec={k:r.get(k,'') for k in ['name','lifestyle','ecology','order','published','source']}
        if res:
            rec.update({'portal':res['portal'],'org_id':res['org_id'],
                        'proteome_fid':res['proteome_fid'],'proteome_name':res['proteome_name'],
                        'proteome_md5':res['proteome_md5'],'proteome_size':res['proteome_size'],
                        'assembly_fid':res['assembly_fid'],'assembly_name':res['assembly_name'],
                        'assembly_md5':res['assembly_md5'],'assembly_size':res['assembly_size'],'status':'OK'})
            ok+=1
        else:
            rec.update({'portal':portal,'status':'UNRESOLVED'}); unresolved.append(name)
        w.writerow(rec); out.flush()
        if i%25==0: print('...%d/%d  ok=%d'%(i,len(rows),ok), flush=True)
    out.close()
    open(BASE+'unresolved.txt','w').write('\n'.join(unresolved))
    print('DONE  total=%d  OK=%d  UNRESOLVED=%d'%(len(rows),ok,len(unresolved)))
    if unresolved: print('unresolved:', ', '.join(unresolved[:40]))
