#!/usr/bin/env bash
set -uo pipefail
PROJ=~/projects/aicomp-h1h3
OF=$PROJ/orthofinder; IN=$OF/input
rm -rf $IN $OF/run; mkdir -p $IN
n=0
while read p; do
  [ -z "$p" ] && continue
  cp $PROJ/clean/$p.faa $IN/$p.fa && n=$((n+1))
done < $PROJ/qc/qc_pass.txt
echo "OrthoFinder input proteomes: $n"
~/miniforge3/bin/mamba run -n orthofinder orthofinder -f $IN -t 60 -a 6 -o $OF/run 2>&1 | tail -25
echo "ORTHOFINDER_DONE"
