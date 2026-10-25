# Baseline comparison: Terrace, Aspen, Spruce, GraphBolt/DZiG vs GIR

Every system runs the same protocol on the same inputs, written by GIR's benchmark:

```bash
cd gir_constructor
./bench_gir --datasets livejournal     --ops 100000 --export ../data/protocol   # directed
./bench_gir --datasets livejournal_sym --ops 50000  --export ../data/protocol   # symmetrized
```

`data/protocol/<dataset>/` then holds the 50% initial graph (`loaded.el`), the full graph
(`full.el`), the insertion / deletion / mixed streams (100,000 directed edge operations each), the
16 BFS/SSSP sources and the 64 BC sources. SSSP weights are GIR's deterministic `edgeWeight(u,v)`
in [1, 100], recomputed by `common/protocol.h`. All runs use 64 threads unless noted.

| Dir | System | Build | Driver measures |
|---|---|---|---|
| `terrace/` | Terrace (SIGMOD'21) | OpenCilk, `make CILK=1` with `-fopencilk` (3-line Makefile edit in the Terrace repo) | batch insert (sorted, deduplicated), parallel `remove_edge` (Terrace has no batch delete), BFS, SSSP (Bellman-Ford), BC |
| `aspen/` | Aspen (PLDI'19) | stock Makefile, homegrown scheduler, `NUM_THREADS` | `insert_edges_batch` / `delete_edges_batch`, BFS, BC. Symmetric graphs only, no SSSP |
| `spruce/` | Spruce (SIGMOD'24) | existing CMake Release build | concurrent `InsertEdge` / `DeleteEdge` per batch (its intended model); BFS / Delta-stepping SSSP / Brandes BC implemented in the driver over `get_neighbours`, since Spruce ships no algorithms |
| `dzig/` | GraphBolt with DZiG runtime (EuroSys'21) | OpenCilk; three missing `return` statements fixed in the GraphBolt repo (the July binaries hung in their initial computation) | DZiG's own streaming mode: per batch, mutation time and incremental BFS / SSSP recompute from one source |
| `gir/` | GIR in DZiG's streaming setting | links `gir_constructor` objects (and jemalloc if installed) | per batch: `applyBatch` time, then BFS / SSSP from scratch from the same source; with a trailing `lazy` argument `ingestBatch` instead (system `GIR-lazy`), plus one `sync()` per stream |

Run: `terrace/`, `aspen/`, `spruce/`, `gir/`, `dzig/prepare` each have a Makefile; `run_all_remaining.sh`
shows the exact command lines used. Results: `experiments/results/baselines.csv` (one row per
measurement: dataset, system, metric, mode, value, unit, count, threads) and per-run logs in
`experiments/results/baselines_logs/`. GIR's own protocol numbers are in
`experiments/results/gir_benchmark.csv` (directed) and `gir_benchmark_sym.csv` (symmetrized).

## Fairness notes

- **Directed vs symmetric.** Terrace's stock BFS/SSSP/BC and all of Aspen's algorithms are only
  correct on symmetric graphs (their dense edge_map reads out-edges as in-edges). On the directed
  graph the Terrace driver uses push-only BFS / Bellman-Ford (correct, verified by matching reach
  counts) and skips BC; Aspen only runs on the symmetrized graph.
- **Algorithms differ.** GIR's BFS is direction-optimizing; Terrace/Aspen use Ligra-style edgeMap;
  the Spruce driver uses top-down BFS (no reverse edges stored). Terrace's SSSP is Bellman-Ford,
  GIR's and the Spruce driver's are Delta-stepping, DZiG's is KickStarter's incremental relaxation.
- **Update semantics.** GIR maintains WCC/SCC/BiCC on every update; the storage systems
  (Terrace, Aspen, Spruce) only update adjacency; DZiG updates adjacency plus one source's BFS/SSSP
  result. Throughput numbers therefore measure different amounts of work.
- **Batch size** is 10,000 operations for every batched system. Spruce processes each batch with a
  parallel loop of single-edge operations. GIR is reported both single-edge (SPU) and batched (FMU).
- **DZiG BFS and SSSP use the KickStarter engine in every mode**, so "GraphBolt" and "DZiG" are the
  same numbers for these algorithms; DZiG's sparsity-aware runtime only affects GraphBolt-engine
  algorithms (PageRank and similar).
