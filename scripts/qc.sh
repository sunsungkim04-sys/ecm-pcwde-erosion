#!/usr/bin/env bash
# Step 2 QC: BUSCO (proteins mode, fungi_odb10) on proteomes + seqkit assembly stats.
set -uo pipefail
PROJ=~/projects/aicomp-h1h3
QC=$PROJ/qc; MB=~/miniforge3/bin/mamba
mkdir -p $QC/busco $QC/tmp $QC/busco_downloads
SEQKIT=~/miniforge3/envs/prokka/bin/seqkit

echo "[1/3] assembly stats (N50/contig)"
$SEQKIT stats -a -T -j 16 $PROJ/genomes/*.fna.gz > $QC/assembly_stats.tsv 2>/dev/null
echo "  -> $QC/assembly_stats.tsv ($(wc -l < $QC/assembly_stats.tsv) lines)"

echo "[2/3] pre-download BUSCO fungi_odb10 lineage"
$MB run -n busco busco --download_path $QC/busco_downloads --download fungi_odb10 >/dev/null 2>&1

echo "[3/3] BUSCO proteins mode (parallel)"
run_busco(){
  faa=$1; portal=$(basename "$faa" .faa.gz)
  ls $QC/busco/$portal/short_summary*.txt >/dev/null 2>&1 && return
  tmp=$QC/tmp/$portal.faa; zcat "$faa" > "$tmp"
  ~/miniforge3/bin/mamba run -n busco busco -i "$tmp" -l fungi_odb10 -m proteins \
     -c 4 -o "$portal" --out_path $QC/busco --offline --download_path $QC/busco_downloads -f >/dev/null 2>&1
  rm -f "$tmp"
}
export -f run_busco; export QC
ls $PROJ/proteomes/*.faa.gz | xargs -P 14 -I{} bash -c 'run_busco "$@"' _ {}

echo "[done] collating BUSCO summaries"
~/miniforge3/envs/qiime2-amplicon-2024.10/bin/python3 - <<'PY'
import glob, re, os, csv
QC=os.path.expanduser('~/projects/aicomp-h1h3/qc')
rows=[]
for f in glob.glob(QC+'/busco/*/short_summary*.txt'):
    portal=f.split('/busco/')[1].split('/')[0]
    txt=open(f).read()
    m=re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)', txt)
    if m: rows.append(dict(portal=portal, complete=float(m.group(1)), single=float(m.group(2)),
                           dup=float(m.group(3)), frag=float(m.group(4)), missing=float(m.group(5)), n=int(m.group(6))))
rows.sort(key=lambda r:r['complete'])
with open(QC+'/busco_summary.tsv','w',newline='') as o:
    w=csv.DictWriter(o,fieldnames=['portal','complete','single','dup','frag','missing','n'],delimiter='\t'); w.writeheader()
    for r in rows: w.writerow(r)
passed=[r for r in rows if r['complete']>=85]
print('BUSCO summaries: %d | complete>=85%%: %d | worst5: %s'%(len(rows),len(passed),[(r['portal'],r['complete']) for r in rows[:5]]))
PY
