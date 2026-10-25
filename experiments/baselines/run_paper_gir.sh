#!/bin/bash
# GIR steps of run_paper.sh only (rerun after a GIR change). Usage: run_paper_gir.sh <dataset> <sweep> [sym=1]
set -u
D=$1; SWEEP=$2; SYM=${3:-1}; MAXT=64
G=/home/cxl6029/projects/GIR; B=$G/experiments/baselines; P=$G/data/protocol; O=$G/experiments/results/paper/$D
CSV=$O/baselines.csv
echo "[$(date +%H:%M:%S)] $D: GIR rerun"
grep -v -E ",GIR(-lazy)?(-inc)?," $CSV > $CSV.tmp && mv $CSV.tmp $CSV   # (drop the old GIR streaming rows)
(cd $G/gir_constructor && OMP_NUM_THREADS=$MAXT ./bench_gir --datasets $D --ops 100000 --batch 10000 --budget 300 --threads $SWEEP --sources 16 --bc-sources 64 --csv $O/gir.csv > $O/gir.log 2>&1)
(cd $G/gir_constructor && OMP_NUM_THREADS=$MAXT ./bench_gir --datasets $D --ops 100000 --batch 10000 --budget 300 --skip-queries --blocks-on-demand --csv $O/gir_bod.csv > $O/gir_bod.log 2>&1)
[ $SYM = 1 ] && (cd $G/gir_constructor && OMP_NUM_THREADS=$MAXT ./bench_gir --datasets ${D}_sym --ops 50000 --batch 10000 --budget 300 --sources 16 --bc-sources 64 --csv $O/gir_sym.csv > $O/gir_sym.log 2>&1)
for m in "" lazy inc lazy-inc; do OMP_NUM_THREADS=$MAXT $B/gir/gir_stream $P/$D $D $CSV 10000 10 $m > $O/gir_stream_${m:-fmu}.log 2>&1; done
grep -E "^FMU.*MEPS" $O/gir.log $O/gir_bod.log
echo "[$(date +%H:%M:%S)] $D: GIR rerun done"
