#!/usr/bin/env bash
cd ~/projects/aicomp-h1h3
PY=~/miniforge3/envs/qiime2-amplicon-2024.10/bin/python3
R=h1_pcwde/refs
echo "### GH7";  $PY scripts/catalytic_check.py GH7  $R/P62694.fasta 229,234     E,E
echo "### GH6";  $PY scripts/catalytic_check.py GH6  $R/P07987.fasta 245         D
echo "### GH28"; $PY scripts/catalytic_check.py GH28 $R/P26214.fasta 180,201,202 D,D,D
echo "### AA9";  $PY scripts/catalytic_check.py AA9  $R/Q7Z9M7.fasta 20,108      H,H
echo "### GH5";  $PY scripts/catalytic_check.py GH5  $R/P07982.fasta 239,350     E,E
echo ALL_CATALYTIC_DONE
