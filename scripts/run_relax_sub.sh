#!/usr/bin/env bash
# Subsampled RELAX for a large family: subsample trimmed alignment -> tree -> label -> HyPhy RELAX.
# Usage: run_relax_sub.sh <FAM> <N_test> <N_ref>
set -uo pipefail
FAM=$1; NT=${2:-80}; NR=${3:-80}
PROJ=/home1/minseo1101/projects/aicomp-h1h3; C=$PROJ/h1_pcwde/codon
PY=~/miniforge3/envs/qiime2-amplicon-2024.10/bin/python3
echo "[1] subsample $FAM ($NT test + $NR ref)"
$PY $PROJ/scripts/subsample_codon.py $FAM $NT $NR
echo "[2] gene tree $FAM"
~/miniforge3/envs/phylogeny/bin/FastTree -nt -gtr -quiet $C/$FAM.sub.codon.fasta > $C/$FAM.sub.tree.nwk 2>/dev/null
echo "[3] label ECM {Test} $FAM"
$PY - "$FAM" <<'PYEOF'
import sys,re
FAM=sys.argv[1]; C="/home1/minseo1101/projects/aicomp-h1h3/h1_pcwde/codon"
tree=open(f"{C}/{FAM}.sub.tree.nwk").read()
ecm=set(l.strip() for l in open(f"{C}/{FAM}.sub.ecm_tips.txt") if l.strip())
def repl(m):
    name=m.group(1); return name+('{Test}' if name in ecm else '')+':'
tree2=re.sub(r'(?<=[(,])([A-Za-z0-9_.\-]+):', repl, tree)
open(f"{C}/{FAM}.sub.labeled.nwk","w").write(tree2)
print(f"{FAM}: labeled {tree2.count('{Test}')} ECM tips")
PYEOF
echo "[4] HyPhy RELAX $FAM"
~/miniforge3/bin/mamba run -n hyphy hyphy relax --alignment $C/$FAM.sub.codon.fasta --tree $C/$FAM.sub.labeled.nwk --test Test --output $C/$FAM.sub.relax.json 2>&1 | grep -iE "relaxation parameter|K\) =|p =|error|Check" | head -20
echo "RELAX_DONE_${FAM}_SUB"
