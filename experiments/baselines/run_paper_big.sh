#!/bin/bash
# (Variant for the billion-edge graphs: no symmetrized GIR, longer timeouts.)
# Every experiment of the paper's evaluation (Sec. 7) on one dataset, all
# systems, one at a time (the machine is shared; nothing runs concurrently).
# Usage: run_paper.sh <dataset> [threads for the sweep, default 1,2,4,8,16,32,64] [max threads, default 64]
# Needs data/protocol/<dataset> and data/protocol/<dataset>_sym (bench_gir --export),
# and data/protocol/<dataset>/dzig (dzig/prepare). Output: experiments/results/paper/<dataset>/.
set -u
D=$1; SWEEP=${2:-1,2,4,8,16,32,64}; MAXT=${3:-64}
G=/home/cxl6029/projects/GIR; B=$G/experiments/baselines; P=$G/data/protocol
O=$G/experiments/results/paper/$D; mkdir -p $O
CSV=$O/baselines.csv
log() { echo "[$(date +%H:%M:%S)] $D: $*"; }
log "start (sweep $SWEEP, max threads $MAXT)"

# ---- GIR: update modes (SPU, FMU, lazy), BFS/SSSP/BC, BFS thread sweep ----
log "GIR directed"
(cd $G/gir_constructor && OMP_NUM_THREADS=$MAXT ./bench_gir --datasets $D --ops 100000 --batch 10000 --budget 300 \
  --threads $SWEEP --sources 16 --bc-sources 64 --csv $O/gir.csv > $O/gir.log 2>&1)
grep -E "^SPU|^FMU|^LAZY.*sync|bfs \(dir|sssp \(delta|bc \(" $O/gir.log | tail -12
log "GIR blocks on demand"
(cd $G/gir_constructor && OMP_NUM_THREADS=$MAXT ./bench_gir --datasets $D --ops 100000 --batch 10000 --budget 300 \
  --skip-queries --blocks-on-demand --csv $O/gir_bod.csv > $O/gir_bod.log 2>&1)
grep -E "^FMU.*MEPS|^LAZY.*sync" $O/gir_bod.log
# (no symmetrized GIR here: the symmetric twitter / friendster graphs exceed GIR's 32-bit CSR offsets)
log "GIR streaming (from scratch, incremental)"
for m in "" lazy inc lazy-inc; do
  OMP_NUM_THREADS=$MAXT $B/gir/gir_stream $P/$D $D $CSV 10000 10 $m > $O/gir_stream_${m:-fmu}.log 2>&1
done
grep -h "endtoend\|mismatch" $O/gir_stream_*.log | head -30

# ---- Terrace: directed (updates, push-only BFS), symmetric (updates, BFS, BC) ----
log "Terrace directed"
CILK_NWORKERS=$MAXT timeout 14400 $B/terrace/terrace_driver $P/$D $D $CSV 10000 > $O/terrace.log 2>&1; echo "  exit $?"
grep "MEPS\|_ms" $O/terrace.log
log "Terrace symmetric"
CILK_NWORKERS=$MAXT timeout 14400 $B/terrace/terrace_driver $P/${D}_sym ${D}_sym $CSV 10000 > $O/terrace_sym.log 2>&1; echo "  exit $?"
grep "MEPS\|_ms" $O/terrace_sym.log

# ---- Aspen (symmetric only): updates, BFS, BC ----
log "Aspen symmetric"
NUM_THREADS=$MAXT timeout 14400 $B/aspen/aspen_driver $P/${D}_sym ${D}_sym $CSV 10000 > $O/aspen_sym.log 2>&1; echo "  exit $?"
grep "MEPS\|_ms" $O/aspen_sym.log

# ---- Spruce: updates, BFS, SSSP, BC (its mixed stream can hang: timeout, then queries alone) ----
log "Spruce directed"
OMP_NUM_THREADS=$MAXT timeout ${SPRUCE_TIMEOUT:-14400} $B/spruce/spruce_driver $P/$D $D $CSV 10000 > $O/spruce.log 2>&1; rc=$?; echo "  exit $rc"
if [ $rc -ne 0 ]; then
  echo "$D,Spruce,mixed_MEPS,TIMEOUT-${SPRUCE_TIMEOUT:-14400}s,-1,MEPS,0,$MAXT" >> $CSV
  OMP_NUM_THREADS=$MAXT QUERIES_ONLY=1 timeout 14400 $B/spruce/spruce_driver $P/$D $D $CSV 10000 > $O/spruce_q.log 2>&1; echo "  queries exit $?"
fi
grep -h "MEPS\|_ms" $O/spruce*.log

# ---- DZiG / GraphBolt (KickStarter engine): BFS and SSSP streams ----
log "DZiG"
timeout 14400 $B/dzig/run_dzig.sh $P/$D $D $CSV $MAXT 10000 10 > $O/dzig.log 2>&1; echo "  exit $?"
tail -24 $O/dzig.log

# ---- BFS thread sweep for the baselines (queries only) ----
for t in ${SWEEP//,/ }; do
  [ $t -eq $MAXT ] && continue  # (the full runs above are the max-thread point)
  log "sweep $t threads"
  CILK_NWORKERS=$t QUERIES_ONLY=1 NO_BC=1 timeout 14400 $B/terrace/terrace_driver $P/$D $D $CSV 10000 > $O/sweep_terrace_$t.log 2>&1
  NUM_THREADS=$t QUERIES_ONLY=1 NO_BC=1 timeout 14400 $B/aspen/aspen_driver $P/${D}_sym ${D}_sym $CSV 10000 > $O/sweep_aspen_$t.log 2>&1
  OMP_NUM_THREADS=$t QUERIES_ONLY=1 NO_BC=1 NO_SSSP=1 timeout 14400 $B/spruce/spruce_driver $P/$D $D $CSV 10000 > $O/sweep_spruce_$t.log 2>&1
done
log "done"
