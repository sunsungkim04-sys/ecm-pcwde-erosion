#!/usr/bin/env python3
# Parse mycorrhizal species list (.md) + Miyauchi DS-2 -> pre-manifest with portal codes where known.
import pandas as pd, re, warnings
warnings.filterwarnings('ignore')
BASE='/home1/minseo1101/projects/aicomp-h1h3/species_lists/'

def norm(s):
    s=str(s).lower()
    s=re.sub(r'\bv[\d.]+\b',' ',s)          # drop version tokens like v1.0/v2.0
    s=re.sub(r'[^a-z0-9 ]',' ',s)
    s=re.sub(r'\s+',' ',s).strip()
    return s

lifemap={'Ectomycorrhizae':'ECM','Arbuscular mycorrhizae':'AM','Ericoid mycorrhizae':'ERM',
         'Orchid mycorrhizae':'ORM','Saprotroph':'SAP','Wood decayer':'WD','Pathogen':'PATH',
         'Endophyte':'ENDO','Yeast':'YEAST','Parasite':'PARA','Lichen symbiont':'LICHEN'}

# --- DS-2 (has portal codes) ---
ds2=pd.read_excel(BASE+'DS2.xlsx', header=None, skiprows=5)[[0,1,4,7,8]]
ds2.columns=['portal','species','order','ecology','strains']
ds2=ds2.dropna(subset=['portal']).reset_index(drop=True)
ds2['norm_strain']=ds2['strains'].map(norm)
ds2['norm_species']=ds2['species'].map(norm)
strain2portal=dict(zip(ds2['norm_strain'],ds2['portal']))
# species-level fallback only when unique
spcount=ds2['norm_species'].value_counts()
species2portal={k:v for k,v in zip(ds2['norm_species'],ds2['portal']) if spcount[k]==1}

# --- mycorrhizal .md table ---
rows=[]
for line in open(BASE+'species_myco.md',encoding='utf-8'):
    line=line.strip()
    if not line.startswith('|'): continue
    cells=[c.strip() for c in line.strip('|').split('|')]
    if len(cells)!=7: continue
    eco,order,name,size,genes,pub,m2020=cells
    if eco=='Ecology' or set(eco)<=set('- '): continue
    rows.append(dict(name=name,ecology=eco,order=order,size=size,genes=genes,published=pub,m2020=m2020))
myco=pd.DataFrame(rows)
myco['norm']=myco['name'].map(norm)

def match_portal(n):
    if n in strain2portal: return strain2portal[n]
    sp=' '.join(n.split()[:2])
    if sp in species2portal: return species2portal[sp]
    return ''
myco['portal']=myco['norm'].map(match_portal)
myco['lifestyle']=myco['ecology'].map(lambda e: lifemap.get(e,'OTHER'))
myco['source']='mycocosm_list'
myco['needs_resolution']=myco['portal']==''

# --- comparison group (non-mycorrhizal from DS-2) ---
comp=ds2[~ds2['ecology'].str.contains('mycorrhiz',case=False,na=False)].copy()
comp['name']=comp['strains']; comp['published']='Miyauchi 2020'; comp['m2020']='included'
comp['lifestyle']=comp['ecology'].map(lambda e: lifemap.get(e,'OTHER'))
comp['source']='ds2_comparison'; comp['needs_resolution']=False

keep=['name','species' if False else 'name','lifestyle','ecology','order','portal','published','source','needs_resolution']
cols=['name','lifestyle','ecology','order','portal','published','source','needs_resolution']
allsp=pd.concat([myco[cols], comp[cols]], ignore_index=True)
allsp.to_csv(BASE+'manifest_pre.tsv', sep='\t', index=False)

print('=== MYCORRHIZAL (from .md) ===', len(myco))
print(myco['lifestyle'].value_counts().to_string())
print('  portal matched:', (myco.portal!='').sum(), '| need JGI resolution:', myco.needs_resolution.sum())
print('=== COMPARISON (from DS-2) ===', len(comp))
print(comp['lifestyle'].value_counts().to_string())
print('=== TOTAL pre-manifest ===', len(allsp), '| need resolution:', allsp.needs_resolution.sum())
print('  saved -> manifest_pre.tsv')
