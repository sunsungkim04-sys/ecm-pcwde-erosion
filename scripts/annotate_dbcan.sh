#!/usr/bin/env bash
set -uo pipefail
PROJ=/home1/minseo1101/projects/aicomp-h1h3
DB=/home1/minseo1101/pk6_analysis/dbcan/db
OUT=$PROJ/annotation/dbcan; mkdir -p $OUT
run_one(){
  faa=$1; p=$(basename "$faa" .faa)
  DB=/home1/minseo1101/pk6_analysis/dbcan/db
  OUT=/home1/minseo1101/projects/aicomp-h1h3/annotation/dbcan
  [ -s $OUT/$p/overview.tsv ] && return
  ~/miniforge3/bin/mamba run -n dbcan run_dbcan CAZyme_annotation \
     --input_raw_data "$faa" --mode protein --output_dir $OUT/$p \
     --db_dir $DB --methods hmm --threads 3 >/dev/null 2>&1
}
export -f run_one
# run on QC-passed proteomes only
mkdir -p $PROJ/logs
awk '{print "/home1/minseo1101/projects/aicomp-h1h3/clean/"$1".faa"}' $PROJ/qc/qc_pass.txt \
  | xargs -P 10 -I{} bash -c 'run_one "$@"' _ {}
echo "DBCAN_DONE  outputs: $(ls -d $OUT/*/ 2>/dev/null | wc -l)"
