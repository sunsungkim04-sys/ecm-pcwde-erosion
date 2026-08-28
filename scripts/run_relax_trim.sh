#!/usr/bin/env bash
# Trimmed RELAX: drop gap-heavy codon columns, then tree + label ECM {Test} + HyPhy RELAX.
# Usage: run_relax_trim.sh <FAM> [maxgap=0.5]
set -uo pipefail
FAM=$1; MAXGAP=${2:-0.5}
PROJ=/home1/minseo1101/projects/aicomp-h1h3; C=$PROJ/h1_pcwde/codon
PY=~/miniforge3/envs/qiime2-amplicon-2024.10/bin/python3
echo "[1] trim codon columns (maxgap=$MAXGAP)  $FAM"
$PY $PROJ/scripts/trim_codon.py $C/$FAM.codon.fasta $C/$FAM.trim.codon.fasta $MAXGAP
echo "[2] gene tree (FastTree -nt) $FAM"
~/miniforge3/envs/phylogeny/bin/FastTree -nt -gtr -quiet $C/$FAM.trim.codon.fasta > $C/$FAM.trim.tree.nwk 2>/dev/null
echo "[3] label ECM tips {Test} $FAM"
$PY $PROJ/scripts/label_tree_trim.py $FAM
echo "[4] HyPhy RELAX $FAM"
~/miniforge3/bin/mamba run -n hyphy hyphy relax --alignment $C/$FAM.trim.codon.fasta --tree $C/$FAM.trim.labeled.nwk --test Test --output $C/$FAM.trim.relax.json 2>&1 | grep -iE "relaxation|k =|p-value|p =|Test |Reference|Likelihood|error|Check" | head -30
echo "RELAX_DONE_$FAM"
