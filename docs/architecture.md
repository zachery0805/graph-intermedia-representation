# Architecture Notes

Quick map of `gir_constructor/` for orienting future changes. This is not exhaustive API docs — see
the code itself and `README.md` for the public `GIR` class API.

## Core system (`gir_constructor/`)

- `gir.h` / `gir.cpp` — the `GIR` class: builds and updates SCC / WCC / Block-Cut-Tree decompositions.
- `graph.h` / `graph.cpp` — directed CSR graph.
- `graph_undirected.h` / `graph_undirected.cpp` — undirected CSR graph.
- `bicc_bfs.h` — biconnected components via BFS.
- `block_cut_tree.h` — Block-Cut Tree structure.
- `dag_graph.h` — SCC-condensed DAG.
- `main.cpp` — entry point. `gir_test.cpp` — test suite (BFS, SSSP, edge updates), compiled into the main binary.

## `includes/` — algorithm library, grouped by family

**K-shortest-paths (KSP) family** — Yen's algorithm and PeeK variants (PeeK = incremental K-shortest-paths):
`KSP.h` (base class) → `Yen.h`, `OptYen.h` (incremental/optimized Yen's), `PeeKAdaptiveWithEdgeSwap.h`,
`PeeKAdaptiveWithEdgeSwapParallel.h`, `PeeKAdaptiveWithStatusArray.h`, `PeeKEdgeSwap.h`,
`PeeKRegeneration.h`, `PeeKStatusArray.h`.

**SSSP / Delta-Stepping family**: `DeltaStepping.h` (templated base) → `DeltaSteppingStatic.h`
(concrete specialization); supporting structures `BucketListStatic.h`, `RequestListStatic.h`,
`Request.h`, `Partition.h`, `SsspTree.h`, `MaxHeapNoDuplicate.h`.

**Graph I/O / representation**: `BasicGraph.h`, `GraphRW.h` (file I/O), `KSPGraph.h` (edge-removal
overlay for KSP), `NodeList.h` (parallel frontier buffer).

**Misc core types**: `Node.h` (NODE_ID/EDGE_ID typedefs, used almost everywhere), `Profiling.h`
(RAII scoped timer), `Statistics.h`.

Note: the KSP/PeeK/Delta-Stepping cluster (attributed to Alex Schickedanz, Uni Frankfurt, 2017/2018 in
file headers) looks like an adopted external K-Shortest-Paths codebase integrated alongside GIR's own
SCC/WCC/BCT logic — treat it as a semi-external library subgroup rather than assuming it follows the
same conventions as the rest of the code.

**SCC/WCC baseline implementations**: `includes/scc_cpu/`, `includes/wcc_cpu/` — batch-testable via
`batch_test*.sh` in each directory (see README "Tests" section).

## Related but separate trees

- `graph_converter/` — format conversion utilities (edge-list ↔ CSR binary ↔ Matrix Market etc.),
  independent of `gir_constructor/`.
- `standard/iCC/` — a separate baseline/reference suite (SCC/WCC/BiCC/bridge detection, CPU+GPU
  variants), used for comparison, not part of the GIR system itself.
- `legacy/` — pre-`gir_constructor` prototypes, superseded, kept only for historical reference.
