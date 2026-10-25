#!/bin/bash
# GraphBolt/DZiG (KickStarter engine for BFS/SSSP) on the GIR protocol streams.
# Usage: run_dzig.sh <protocol dir> <dataset label> <csv> [threads=64] [batch=10000] [batches=10]
# Needs: prepare <protocol dir> <protocol dir>/dzig  (adjacency + stream files).
set -e
P=$1; DS=$2; CSV=$3; TH=${4:-64}; B=${5:-10000}; NB=${6:-10}
GB=/home/cxl6029/projects/GIR/experiments/dzig_bench/graphbolt
MI=$GB/lib/mimalloc/out/release
export LD_LIBRARY_PATH=$MI LD_PRELOAD=$MI/libmimalloc.so CILK_NWORKERS=$TH
SRC=$(head -1 $P/sources.txt)
LOG=/home/cxl6029/projects/GIR/experiments/results/baselines_logs
mkdir -p $LOG
for alg in BFS SSSP_w; do
  for s in insert delete mixed; do
    if [ $alg = SSSP_w ]; then adj=$P/dzig/loaded_w.adj; st=$P/dzig/${s}_w.stream; else adj=$P/dzig/loaded.adj; st=$P/dzig/${s}.stream; fi
    out=$LOG/dzig_${DS}_${alg}_${s}.log
    $GB/apps/$alg -source $SRC -fixedBatchSize -nEdges $B -numberOfUpdateBatches $NB -streamPath $st $adj > $out 2>&1
    python3 - "$out" "$CSV" "$DS" "$alg" "$s" "$TH" "$B" <<'PY'
import sys, re
log, csv, ds, alg, stream, th, B = sys.argv[1:]
txt = open(log).read()
init = float(re.search(r"Initial graph processing : ([0-9.e-]+)", txt).group(1))
dels = [float(x) for x in re.findall(r"Edge deletion time : ([0-9.e-]+)", txt)]
adds = [float(x) for x in re.findall(r"Edge addition time : ([0-9.e-]+)", txt)]
fin = [float(x) for x in re.findall(r"Finished batch : ([0-9.e-]+)", txt)]
nb = len(fin); ops = nb * int(B)
algname = "bfs" if alg == "BFS" else "sssp"
rows = [
  (f"{algname}_initial_ms", "from-scratch", init * 1000, "ms", 1),
  (f"{stream}_mutation_MEPS", f"{algname}_batch_{B}", ops / 1e6 / (sum(dels) + sum(adds)), "MEPS", ops),
  (f"{algname}_{stream}_recompute_ms", f"incremental_batch_{B}", 1000 * sum(fin) / nb, "ms", nb),
  (f"{algname}_{stream}_endtoend_ms", f"batch_{B}", 1000 * (sum(dels) + sum(adds) + sum(fin)) / nb, "ms", nb),
]
with open(csv, "a") as f:
    for m, mode, v, u, c in rows:
        f.write(f"{ds},DZiG,{m},{mode},{v:.6g},{u},{c},{th}\n")
        print(f"  {m:28s} {mode:24s} {v:.6g} {u}")
PY
  done
done
