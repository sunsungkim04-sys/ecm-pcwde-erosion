#!/usr/bin/env bash
# Build gene tree + label ECM + run HyPhy RELAX for a PCWDE family. Usage: run_relax.sh GH7
set -uo pipefail
FAM=$1
PROJ=/home1/minseo1101/projects/aicomp-h1h3; C=$PROJ/h1_pcwde/codon
PY=~/miniforge3/envs/qiime2-amplicon-2024.10/bin/python3
echo "[1] codon alignment"; $PY $PROJ/scripts/codon_align.py $FAM
echo "[2] gene tree (FastTree -nt)"; ~/miniforge3/envs/phylogeny/bin/FastTree -nt -gtr -quiet $C/$FAM.codon.fasta > $C/$FAM.tree.nwk 2>/dev/null
echo "[3] label ECM tips {Test}"; $PY $PROJ/scripts/label_tree.py $FAM
echo "[4] HyPhy RELAX"; ~/miniforge3/bin/mamba run -n hyphy hyphy relax --alignment $C/$FAM.codon.fasta --tree $C/$FAM.labeled.nwk --test Test --output $C/$FAM.relax.json 2>&1 | grep -iE "relaxation|k =|p-value|p =|Test |Reference|Likelihood" | head -25
echo "RELAX_DONE_$FAM"
