# GIR — Graph Incremental Reachability

GIR is a C++ system for maintaining and querying graph structures (Strongly Connected Components, Weakly Connected Components, and Block-Cut Trees) on large-scale graphs that undergo dynamic edge insertions and deletions. It targets real-world graphs with hundreds of millions of edges and uses OpenMP for parallelism.

## Features

- **Dynamic graph structures** — maintain SCC/WCC decompositions and Block-Cut Trees under edge insertions and deletions
- **Parallel algorithms** — OpenMP-based BFS, Delta-Stepping SSSP, and component detection
- **K-Shortest Paths** — Yen's algorithm and several PeeK-Adaptive variants
- **Multiple graph formats** — CSR binary and edge-list text; conversion utilities included
- **Benchmark datasets** — LiveJournal, Twitter, RoadUSA, web-Google, YouTube, and more

## Repository Layout

```
GIR/
├── gir_constructor/       # Core system (build and run this)
│   ├── gir.h / gir.cpp    # GIR class — builds and updates graph decompositions
│   ├── graph.h / .cpp     # Directed CSR graph
│   ├── graph_undirected.* # Undirected CSR graph
│   ├── bicc_bfs.h         # Biconnected components via BFS
│   ├── block_cut_tree.h   # Block-Cut Tree structure
│   ├── dag_graph.h        # SCC-condensed DAG
│   ├── main.cpp           # Entry point
│   ├── gir_test.cpp       # Test suite (BFS, SSSP, edge updates)
│   ├── Makefile
│   ├── build_and_run.sh
│   └── includes/          # Algorithm libraries (SCC, WCC, SSSP, KSP)
│       ├── scc_cpu/
│       ├── wcc_cpu/
│       ├── DeltaStepping.h
│       ├── KSP.h / Yen.h
│       └── PeeKAdaptive*.h
├── graph_converter/       # Graph format conversion utilities (C++ and Python)
│   ├── edge_list_to_csr_bin/
│   └── *.py
├── standard/              # Baseline / reference implementations (iCC suite), for comparison only
├── legacy/                # Superseded pre-gir_constructor prototypes, kept for history only
├── experiments/           # Experiment results (results/) and analysis (analysis/) — see docs/experiment-plan.md
├── docs/                  # Datasets manifest, experiment plan, architecture notes
└── data/                  # Benchmark graph datasets (binary CSR format) — not tracked in git except example_graph/
    ├── example_graph/     # Small 12-vertex toy graph
    ├── livejournal/       # ~4.8M vertices, ~69M edges
    ├── twitter/
    ├── RoadUSA_sym/
    └── ...
```

See [`docs/architecture.md`](docs/architecture.md) for a map of `gir_constructor/includes/`,
[`docs/datasets.md`](docs/datasets.md) for the full dataset list, and
[`docs/experiment-plan.md`](docs/experiment-plan.md) for what's planned experimentally.

## Prerequisites

| Requirement | Version |
|---|---|
| C++ compiler | g++ with C++17 support |
| OpenMP | Any standard version |
| Make | GNU make |
| Python | 3.x (graph conversion utilities only) |

## Build

```bash
cd gir_constructor
make clean && make all
```

Optional build flags:

```bash
make DEBUG=1        # enable debug output
make OUTPUT_TIME=1  # print timing information
```

Or use the convenience script:

```bash
./build_and_run.sh --debug --time
```

## Preparing Graph Data

GIR reads graphs in **CSR binary format** — four files per graph:

| File | Description |
|---|---|
| `fw_begin.bin` | Forward adjacency row pointers |
| `fw_adjacent.bin` | Forward adjacency column indices |
| `bw_begin.bin` | Backward (transpose) row pointers |
| `bw_adjacent.bin` | Backward adjacency column indices |

### Convert an edge-list text file to CSR binary

```bash
cd graph_converter/edge_list_to_csr_bin
make
./converter <input_edge_list.txt> <output_dir/>
```

### Other conversion utilities (Python)

```bash
# Edge list → Matrix Market
python3 graph_converter/edge_list_to_matrix_market.py input.txt output.mtx

# See graph_converter/ for the full set of conversion scripts
```

## Running

Edit the data path in `gir_constructor/main.cpp` before building:

```cpp
std::string path_undirected = "/your/path/to/data/livejournal/";
std::string path_directed   = "/your/path/to/data/livejournal/";
```

Then run:

```bash
./gir
```

A small 12-vertex example graph is included in `data/example_graph/` for quick testing.

## GIR API

```cpp
#include "gir.h"

GIR gir(path_directed, path_undirected);

// Build decompositions
gir.buildWCC();           // Weakly Connected Components
gir.buildSCC_DAG();       // Strongly Connected Components + condensed DAG
gir.buildBCT();           // Block-Cut Tree (biconnected components)

// Query
auto& sccs  = gir.getSCCs();
auto& wccs  = gir.getWCCs();
auto& dag   = gir.getDAG();
int scc_id  = gir.getNodeToScc(v);
int wcc_id  = gir.getNodeToWcc(v);

// Dynamic updates
gir.updateEdgeDirected(u, v, /*insert=*/true);
gir.updateEdgeDirected(u, v, /*insert=*/false);
gir.updateEdgeUndirected(u, v, insert);
```

## Tests

The test suite in `gir_test.cpp` is compiled into the main binary:

| Test | Description |
|---|---|
| `ParallelBFS` | Multi-source BFS on the directed graph |
| `BCT_BFS_Test` | BFS on the Block-Cut Tree |
| `EdgeUpdateTest` | Batch edge insertions and deletions |
| `ParallelSSSP` | Delta-Stepping single-source shortest paths |

SCC and WCC algorithm batch tests are in `includes/scc_cpu/batch_test*.sh` and `includes/wcc_cpu/batch_test*.sh`.

## Benchmark Datasets

| Dataset | Vertices | Edges | Notes |
|---|---|---|---|
| example_graph | 12 | — | Toy graph for quick tests |
| com-youtube | ~1.1M | ~3M | YouTube social network |
| web-Google | ~875K | ~5M | Web crawl graph |
| soc-LiveJournal1 | ~4.8M | ~69M | LiveJournal social network |
| RoadUSA | — | — | US road network |
| Twitter | — | — | Large-scale Twitter follow graph |

Datasets are stored as CSR binary files under `data/`.

## License

See `LICENSE` if present, or contact the repository owner.
