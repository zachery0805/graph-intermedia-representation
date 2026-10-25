#!/bin/bash
# Runs, one at a time: DZiG + GIR streaming (directed), then the symmetric
# comparison (Aspen, Terrace, GIR). Results append to experiments/results/baselines.csv.
G=/home/cxl6029/projects/GIR
B=$G/experiments/baselines; R=$G/experiments/results/baselines.csv; P=$G/data/protocol; L=$G/experiments/results/baselines_logs
echo "== DZiG"; $B/dzig/run_dzig.sh $P/livejournal livejournal $R 64 10000 10
echo "== GIR stream"; OMP_NUM_THREADS=64 $B/gir/gir_stream $P/livejournal livejournal $R 10000 10 > $L/gir_stream_lj.log 2>&1; grep MEPS\|_ms $L/gir_stream_lj.log
echo "== Aspen sym"; NUM_THREADS=64 $B/aspen/aspen_driver $P/livejournal_sym livejournal_sym $R 10000 > $L/aspen_ljsym.log 2>&1; grep "MEPS\|_ms\|reached" $L/aspen_ljsym.log
echo "== Terrace sym"; CILK_NWORKERS=64 $B/terrace/terrace_driver $P/livejournal_sym livejournal_sym $R 10000 > $L/terrace_ljsym.log 2>&1; grep "MEPS\|_ms\|reached" $L/terrace_ljsym.log
echo "== GIR sym"; cd $G/gir_constructor && OMP_NUM_THREADS=64 ./bench_gir --datasets livejournal_sym --ops 50000 --batch 10000 --budget 300 --sources 16 --bc-sources 64 --csv $G/experiments/results/gir_benchmark_sym.csv > $L/gir_ljsym.log 2>&1; grep -E "SPU:|FMU|bfs|sssp|bc " $L/gir_ljsym.log
echo ALL DONE
