#!/usr/bin/env python3
# Retry UNRESOLVED species: DS-2 portal lookup / strain query / portal-code-from-filenames fallback.
import sys, csv, re, collections
sys.path.insert(0,'/home1/minseo1101/projects/aicomp-h1h3/scripts')
import resolve_portals as R
import pandas as pd, warnings; warnings.filterwarnings('ignore')
BASE='/home1/minseo1101/projects/aicomp-h1h3/species_lists/'

# DS-2 binomial -> portal
ds2=pd.read_excel(BASE+'DS2.xlsx', header=None, skiprows=5)[[0,1]]
ds2.columns=['portal','species']; ds2=ds2.dropna(subset=['portal'])
binom2portal={}
for _,r in ds2.iterrows():
    binom2portal.setdefault(R.norm(r['species']), r['portal'])

def extract_portal(d):
    c=collections.Counter()
    for org in d.get('organisms',[]):
        for f in org.get('files',[]):
            m=re.match(r'^([A-Z][A-Za-z]{2,}\d+)_', f.get('file_name',''))
            if m: c[m.group(1)]+=1
    return [p for p,_ in c.most_common(3)]

def resolve2(name):
    n=R.norm(name); binom=' '.join(n.split()[:2])
    portal=binom2portal.get(binom)
    tries=[]
    if portal: tries.append(portal)
    tries += [n, binom]
    for q in tries:
        res=R.resolve_query(q, want_portal=portal, want_name=name)
        if res: return res
    # fallback: portal codes from filenames
    for sq in (n, binom, name):
        d=R.search(sq)
        if not d: continue
        for p in extract_portal(d):
            res=R.resolve_query(p, want_portal=p, want_name=name)
            if res: return res
    return None

rows=list(csv.DictReader(open(BASE+'manifest.tsv'), delimiter='\t'))
fields=rows[0].keys()
fixed=0
for r in rows:
    if r.get('status')=='OK': continue
    res=resolve2(r['name'])
    if res:
        r.update({'portal':res['portal'],'org_id':res['org_id'],
                  'proteome_fid':res['proteome_fid'],'proteome_name':res['proteome_name'],
                  'proteome_md5':res['proteome_md5'],'proteome_size':res['proteome_size'],
                  'assembly_fid':res['assembly_fid'],'assembly_name':res['assembly_name'],
                  'assembly_md5':res['assembly_md5'],'assembly_size':res['assembly_size'],'status':'OK'})
        fixed+=1; print('FIXED %-34s -> %s'%(r['name'][:34], res['portal']), flush=True)
    else:
        print('STILL %-34s'%r['name'][:34], flush=True)
w=csv.DictWriter(open(BASE+'manifest.tsv','w',newline=''), fieldnames=list(fields), delimiter='\t')
w.writeheader(); [w.writerow(r) for r in rows]
still=[r['name'] for r in rows if r.get('status')!='OK']
open(BASE+'unresolved.txt','w').write('\n'.join(still))
print('RETRY DONE  fixed=%d  still_unresolved=%d'%(fixed,len(still)))
