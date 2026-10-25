#!/bin/bash
G=/home/cxl6029/projects/GIR; R=$G/experiments/results/baselines.csv; P=$G/data/protocol; L=$G/experiments/results/baselines_logs
cd $G/experiments/baselines/spruce
for t in 64 32; do SPRUCE_SKIP_MIXED=1 OMP_NUM_THREADS=$t ./spruce_driver $P/livejournal livejournal $R 10000 > $L/spruce_lj_$t.log 2>&1; done
echo SPRUCE DONE >> $L/spruce_lj_32.log
