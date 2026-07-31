# Changelog

Milestone-level summary of major changes. Day-to-day detail lives in `git log`; this file is for
skimming "what happened when" while writing up the paper.

## 2026-07-29 — Fixed `bfsTaskParallel` being consistently slower than sequential BFS

The first-version benchmark (2026-07-27 entry below) showed `bfsTaskParallel` 2-16.6x *slower* than a
plain sequential BFS on every real dataset tested, worst on the largest ones (twitter: 16.6x slower;
RoadUSA: 11.7x; grid: 15.9x). Root-caused and fixed four distinct, compounding issues, all in
`GIR::bfsTaskParallel` (`gir.cpp`), validated via the existing differential test
(`GirBfsTest::runDifferentialBfsTest`, matches sequential BFS exactly on every dataset after every fix):

1. **`nodeToScc`/`reachable_sccs` lookups on the hot path.** Every edge examined during the BFS called
   `nodeToScc.find(w)` (a `std::map`, i.e. a red-black tree -- O(log V), pointer-chasing, cache-hostile)
   and `reachable_sccs.count(...)` (an `unordered_set` hash lookup). Fixed by snapshotting both into flat
   arrays (`scc_of[v]`, `reachable_sccs[scc_id]` as a `vector<char>`) for O(1) access -- see point 4 below
   for how that snapshot itself had to be scoped correctly.
2. **Block-grouping degenerates to one task per vertex on real graphs.** The implementation partitioned
   each level's frontier by `(scc, block)` and dispatched one OpenMP task per group, mirroring the
   paper's "task parallelism across blocks" model. But most SCCs on real graphs are near-singletons with
   no multi-vertex BiCC block (971K SCCs on livejournal's 4.8M vertices, 8M on twitter's 41.6M) --
   vertices with no block info fall into their own singleton group, so grouping effectively became "one
   OpenMP task per vertex". Diagnostic timing (temporary, since removed) showed the "parallel" region on
   livejournal taking as long as the *entire* sequential BFS baseline despite 128 cores being available --
   dynamic-scheduling dispatch overhead per task exceeded the ~14-edges-per-vertex of real work each task
   did. **Decision, discussed with and confirmed by the user:** dropped block-based grouping entirely in
   favor of plain data parallelism directly over the frontier array (`schedule(dynamic, 256)`). This is a
   real scope change from the paper's stated task+data two-level model -- the implementation no longer
   does "task parallelism across blocks", just data parallelism across vertices. Documented here rather
   than silently: on graphs this fragmented, block-level task granularity isn't a win, and a paper
   claiming this design should either qualify that or use a size-adaptive grouping strategy, which this
   codebase does not attempt.
3. **Fixed OpenMP parallel-region overhead on small frontiers.** Even after fix 2, unconditionally
   parallelizing every level was a net loss on graphs where the reachable region itself is tiny --
   web-Google has 7.77M WCCs out of 8.64M vertices (almost all singleton/near-singleton), so a random BFS
   source typically reaches only a handful of vertices per level, and thread-team spin-up cost alone made
   task-parallel BFS ~19-21x slower than sequential there. Added `if(n >= 2000)` to the `#pragma omp
   parallel for`, which serializes the loop in the calling thread with no thread-team overhead below that
   frontier-size threshold (OpenMP's own `if` clause semantics, not a hand-rolled branch -- no duplicated
   loop body).
4. **The point-1 snapshot cost has to be scoped, not global.** First attempt built `scc_of` by iterating
   *all* of `nodeToScc` (O(V) via `std::map` in-order traversal) once per call. That's fine when the BFS
   is going to visit most of V anyway (twitter/livejournal/RoadUSA: dominated by one giant WCC), but on
   web-Google it made things *worse* than before this whole fix pass (1.28ms -> 283ms, ~220x) since a
   single-vertex BFS was still paying a full 8.64M-entry map traversal. Fixed by computing Stage 1 (SCC-DAG
   reachability) *first* -- its cost is bounded by the reachable region, not V -- then building `scc_of`
   only for vertices in reachable SCCs, using the `sccs[scc_id]` vertex-lists GIR already maintains
   (contiguous vectors, not a map traversal). This bounds every per-call precompute by the actual
   reachable region's size in all cases, and turned out faster even for the giant-WCC datasets (contiguous
   vector iteration beats `std::map` traversal at the same asymptotic cost).

**Measured results** (task-parallel avg / sequential avg / speedup, 20 random sources, `bench_gir`
methodology unless noted):

| Dataset | Before | After | Improvement |
|---|---|---|---|
| livejournal | 17,804ms / 1,235ms / 0.069x (14.4x slower) | 301ms / 1,036ms / **3.44x** (faster than sequential) | 59x |
| RoadUSA | 37,371ms / 3,205ms / 0.086x (11.7x slower) | 6,101ms / 3,421ms / 0.56x (1.78x slower) | 6.1x |
| grid_1000_100000 | 1,955ms / 123ms / 0.063x (15.9x slower) | 323ms / 235ms / 0.73x (1.37x slower) | 6.1x |
| web-Google | 1.28ms / 0.62ms / 0.48x (2.1x slower, original) | 34.1ms / 13.7ms / 0.40x (2.5x slower) | -- (see caveat) |
| twitter | 429,210ms / 25,850ms / 0.060x (16.6x slower) | 8,286ms / 26,392ms / **3.18x** (faster than sequential) | 51.8x |

Caveats, stated plainly rather than glossed over:
- web-Google's "after" number is *worse* than its very first (pre-any-fix) measurement, because that
  original run predates fixes 1-2 entirely (no map/grouping overhead existed yet to fix) -- it's included
  to show the fix-4 regression-and-recovery within this pass, not as a regression from a good baseline.
  web-Google's BFS visits are inherently tiny (fragmented WCC structure), so fixed per-call overhead
  (two O(V)-sized array allocations for `levels`/`scc_of`, `omp_get_max_threads()`, Stage-1 setup) will
  always cost more, in relative terms, than on a dense giant-component graph -- there's a practical floor
  here that further micro-optimization could still chip at but won't eliminate.
- All "after" numbers were measured with a throwaway harness (`quick_bfs_check.cpp`, since removed) that
  runs BFS immediately after construction, before any `insertEdge`/`deleteEdge` calls -- i.e. with the
  dynamic-update overlay (`extra_out_adj`/`removed_directed`) empty. `bench_main.cpp`'s actual benchmark
  order runs BFS *after* the update-throughput phase, so the overlay won't be empty there in practice.
  Fixes 1-4 all still apply regardless of overlay state; only one minor bonus (a raw-CSR fast path that
  bypasses `directedOutNeighbors`'s per-call heap allocation, gated on the overlay being empty) won't
  trigger post-update. The full end-to-end re-benchmark (needed to get final post-update numbers into
  `experiments/results/gir_benchmark.csv`) hasn't been re-run yet -- twitter's update-throughput phase
  alone takes 4+ hours given the still-open Case-2 rank-interval cost (see the 2026-07-27 entries above),
  so this is a real time investment, not a quick re-run.

## 2026-07-27 — Disk cleanup; `data/twitter/` is actually mis-converted (vertex ID space not compacted)

The shared server's disk hit 100% (16GB free out of 7.0T) while investigating adding Friendster
(DZiG's paper cites it at ~2.5B edges / 68.3M vertices) as a benchmark dataset alongside Twitter, to
mirror DZiG's own evaluation methodology.

- **This project's own footprint is not the bottleneck.** `data/` totaled ~153GB out of the 6.6TB used
  disk-wide; the rest is other users'/other projects' usage on this shared machine, which cleanup here
  can't address. Deleted a handful of genuinely redundant/unrelated items to reclaim what's actually
  reclaimable from this project: `data/soc-LiveJournal/` (1.3G, a duplicate/older LiveJournal conversion
  — `data/livejournal/` is kept as canonical), `data/CHEM_5_sym.bin` (433M), `data/GeoLifeNoScale_2.bin`
  (1.1G), `data/com-youtube_sym.bin` (89M, unused by the current benchmark suite), and
  `data/twitter/twitter.txt` (25G raw edge-list text — believed at the time to be safe to drop since the
  converted `.bin` CSR files were assumed correct and complete). Net: 16GB → 43GB free. Adding a
  Friendster-scale dataset (~150-250GB estimated) is on hold pending a real storage solution
  (quota/external storage), not something reclaimable from this project's own data.
- **Correction — `data/twitter/` does NOT match DZiG's Twitter dataset; the `.bin` conversion is broken.**
  Running it through `bench_main.cpp --twitter-only` surfaced `vert_count = 1,468,365,183`, essentially
  equal to `edge_count = 1,468,365,182` (off by exactly 1). A real Twitter follower graph has average
  out-degree in the thousands, so vertex count should be far smaller than edge count (DZiG's Table 2:
  ~41.7M vertices / ~1.5B edges for this exact source, Kwak et al. 2010) — vert≈edge here means the
  vertex ID space was never compacted down from Twitter's large/sparse raw user IDs to a dense `[0, V)`
  range during the original txt→bin conversion. `graph.cpp`'s loader (`vert_count =
  fsize(fw_begin.bin)/sizeof(int) - 1`) is straightforward and not at fault — the `.bin` files themselves
  are wrong. Downstream this produced 1.436B near-singleton "SCCs" (SCC count ≈ vertex count, consistent
  with a corrupted/uncompacted structure) and a `std::bad_array_new_length` crash in `buildBCT()` once
  per-SCC data structures tried to size themselves against 1.4B SCCs. This reopens (rather than closes)
  the scale-mismatch question flagged in the 2026-07-22 entry below — and is now harder to fix in place
  since `twitter.txt` (the raw edge list needed to redo the conversion correctly) was just deleted above;
  no equivalent properly-compacted copy was found elsewhere on this machine (checked `iCC`'s reference
  scripts, which point at a `/mnt/raid0_huge/yuede/...` path that doesn't exist on this host). Twitter is
  dropped from the benchmark's scope until this is resolved.
- **Resolved — a correct, already-compacted copy was sitting unused in the same directory.**
  `data/twitter/twitter.bin` (6.2G, never read by GIR's loader, initially assumed to be another dead
  artifact) turned out to be a self-describing CSR: a 24-byte header of three `int64`s
  (`vert_count=41652231, edge_count=1468365182, total_file_bytes=6206678608`) followed by a `(V+1)`-entry
  `int64` begin-offsets array and an `E`-entry `int32` csr array. Verified internally consistent
  (offsets monotonic non-decreasing, last offset == E, csr values in `[0, V)`, header's declared total
  size matches the actual file size to the byte) and it matches DZiG's Table 2 numbers for Twitter
  exactly (41.7M vertices / 1.5B edges, Kwak et al. 2010) — this is the dataset GIR's loader was supposed
  to have been built from. Wrote `gir_constructor/rebuild_twitter.py` (numpy) to parse this file and
  regenerate the 6-file layout `graph.cpp`/`graph_undirected.cpp` actually read: `fw_begin.bin`/
  `fw_adjacent.bin` as a direct int32 reinterpretation of the header's arrays, `undirected_begin.bin`/
  `undirected_adjacent.bin` as a copy of the same (matches this codebase's existing convention of relying
  on `buildBCT()` to symmetrize an asymmetric undirected CSR rather than shipping a separate symmetrized
  file), and `bw_begin.bin`/`bw_adjacent.bin` built as a proper transpose (sort-edges-by-target, using
  numpy's vectorized `argsort`/`bincount` rather than a scalar scatter loop, since a per-edge Python loop
  over 1.47B edges would've taken hours). Also deleted six further files never read by GIR's loader
  (`head.bin`, `in_degree.bin`, `out_degree.bin`, `bw_head.bin`, `undirected_head.bin`,
  `undirected_out_degree.bin` — ~55G, artifacts of whatever other pipeline produced the broken files).
  Re-running `bench_gir --twitter-only` after the swap now reports `41652231 verts, 1468365182 edges` on
  load — matching the paper. Net effect: `data/twitter/` shrank from 129G to ~24G, and the disk went from
  16GB free to 120GB free over the course of this cleanup+rebuild, entirely from this project's own
  directory (no other users'/projects' data touched). Full experimental numbers for Twitter once the
  benchmark run completes are in `experiments/results/gir_benchmark.csv`.

## 2026-07-22 — Implemented GIR's two-level parallel BFS + first-version benchmarks

`paper/task_parallel.tex` describes GIR's second major contribution (beyond the SPU maintenance
algorithm): a two-level parallel execution model (task parallelism across BC-tree blocks, data
parallelism within each block), with BFS as the worked example. None of this existed in code before —
only construction and dynamic maintenance did. Added `GIR::bfsTaskParallel(source)`:

- **Stage 1** (matches the paper exactly): prune to SCCs reachable from the source's SCC via the
  condensation DAG.
- **Stages 2/3**: the paper's own description ("each task runs BFS from its entry AP, computing local
  levels for all v ∈ B_i") implies a single well-defined entry point per block. That holds *within* one
  SCC's BC-tree (adjacent blocks share exactly one AP, by construction), but **not across SCCs**: two
  different external SCCs can each have their own edge landing on a *different* vertex of the same
  target block, giving that block two genuinely independent entries with no BC-tree relationship
  between them. A differential test against a plain sequential BFS (see below) caught exactly this on
  the 18-vertex `example_graph` — running one local BFS per block from a single "first discovered"
  entry point silently overwrote a shared vertex's correct level with a wrong, larger one. The
  implementation instead runs standard level-synchronous BFS (identical to sequential BFS by
  construction, so this class of bug can't recur) and only uses block membership to *group and dispatch*
  each level's frontier across OpenMP threads — task parallelism across blocks that happen to have
  frontier work in the same level, rather than trying to skip multiple levels inside one block in a
  single shot.
- Added `GirBfsTest` in `gir_test.cpp`: differential correctness check (random sources vs. sequential
  BFS) plus a timing harness.

**First-version experimental data** (GIR only, no external baselines yet — see below for why): construction
time, dynamic-update throughput (MEPS), and BFS runtime (task-parallel vs. sequential) across the real
datasets available under `data/`. Full numbers in `experiments/results/gir_benchmark.csv` and
`docs/experiment-plan.md`.

**Deliberately out of scope for this pass** (flagged, not silently skipped): comparison against the five
external systems the paper's draft tables currently cite (Terrace, Aspen, Spruce, GraphBolt, DZiG) --
three of the five are cloned locally (`~/projects/{terrace,aspen,Spruce}`) but not yet verified to build,
and GraphBolt/DZiG aren't present at all; the paper's two largest datasets (TW/FS at 1--2 SNAP-scale
billion edges) also aren't confirmed to match what's under `data/twitter/` (see the experiment writeup).
SSSP and BC task-parallel engines are likewise not yet implemented (BFS only, per the agreed scope).

## 2026-07-21 — Repository reorganization

- Unified three previously-separate git histories (`gir_constructor/`, `graph_converter/`, plus the
  untracked rest of the tree) into one repo at the top level, preserving all prior commit history.
- Removed a GitHub PAT that had been sitting in plaintext in `gir_constructor`'s git remote URL — no
  longer present in any repo config; the token itself should still be revoked/rotated on GitHub.
- Archived superseded pre-`gir_constructor` prototype files into `legacy/`.
- Added `experiments/` (results + analysis, separate from source) and `docs/` (dataset manifest,
  experiment plan, architecture notes) to support ongoing paper writing.
- Removed regenerable build artifacts and binary test outputs (~6GB in `graph_converter/`, assorted
  `.o`/compiled binaries throughout, a redundant 281MB zip backup of `standard/iCC/`) that had
  accumulated in the working tree.
- Added `.gitignore` covering build artifacts and the `data/` datasets directory (kept on disk, not
  versioned, except the small `example_graph/` toy fixture).

## 2026-07-21 — Deep correctness code review

Full findings, severity, and suggested fixes in [`docs/code-review.md`](docs/code-review.md). Five review
passes across `gir_constructor/` (core GIR class, KSP/Yen/PeeK family, Delta-Stepping/SSSP + graph I/O,
scc_cpu/wcc_cpu baselines) and `graph_converter/`, plus a targeted adversarial re-verification pass on the
highest-impact claims (read directly against source, not just taken on faith from the first pass).

Headline results — **do not treat these as validated until fixed**:
- The dynamic edge-update API (`updateEdgeDirected`/`updateEdgeUndirected` and everything under it in
  `gir.cpp`) is mostly non-functional: several handlers are no-ops or traverse the wrong structure
  entirely. The only test that would have caught this is commented out in `main.cpp`.
- The SCC-condensed DAG is never actually populated — `dag.bin` always records 0 neighbors per SCC, and
  `GIR::getDAG()` always returns empty, regardless of the input graph.
- `graph_converter/edge_list_to_csr_bin/txt_to_bin_undirected.cpp` never symmetrizes edges (the
  reverse-edge insertion is commented out) despite writing to the exact filenames `gir_constructor` loads
  as "undirected" — every "undirected" graph produced by this tool is silently still directed.
  `txt_to_bin_fw.cpp` is a byte-for-byte duplicate of `txt_to_bin_undirected.cpp` and writes to the wrong
  output files entirely, so it never produces a forward CSR at all.
- `scc_cpu`'s SCC-labeling scheme has an ID-space collision (trim-assigned IDs, a hardcoded "largest SCC"
  sentinel, and mice-pivot vertex-ID-based labels all share the same integer range) that can silently
  merge unrelated SCCs.
- `wcc_cpu/wcc_main.cpp` calls `wcc_detection()` with the wrong number/type of arguments and cannot
  currently compile — any assumption that its batch-test results reflect the current `wcc_core.cpp` should
  be treated as unverified.

See `docs/code-review.md` for the full list (including lower-severity dead code, leaks, data races, and a
prioritized fix order) and per-finding suggested fixes. This pass was review-only — no algorithm code was
changed.

## 2026-07-22 — Rewrote GIR's dynamic edge-update (SPU) to match the paper's design

`paper/` (a PPoPP'27 submission draft, not previously in the repo) specifies a five-case
Scope-Prune-Update classification for single-edge insertions/deletions (`paper/maintenance.tex`,
`paper/single_edge_update.tex`) with per-region complexity bounds. The previous `gir.cpp` update
functions (flagged as mostly non-functional in the code review above) have been replaced with an
implementation of that design:

- **WCC layer**: real union-find (`wcc_parent`/`wcc_size`), O(α(n)) merge for Case 1. Case 5c (WCC
  split) uses a lockstep bidirectional BFS to find whichever resulting side is smaller, then rebuilds
  union-find state for the whole affected WCC (bounded by that WCC's size, not the graph's — see the
  "known limitations" note in `gir.h`, this is not yet bounded by the *smaller* side specifically).
- **SCC / condensation-DAG layer**: added a per-SCC topological `rank` (sparse spacing, `rank_to_scc`
  sorted index) so Case 2 can bound cycle detection to the rank interval between the two endpoints
  (bidirectional search + local re-topo-sort of just that interval, not the whole DAG). Case 5b
  (intra-SCC bridge deletion) re-checks strong connectivity with a scoped Tarjan run only over the old
  SCC's own vertices.
- **BiCC/BC-tree layer**: now computed **per SCC** on the induced undirected subgraph, matching the
  paper's formal model in `task_parallel.tex` (`buildBCT()` previously computed one global BCC/BC-tree
  for the whole graph, which silently violated the paper's own definition — fixed as part of this work).
  Cases 3/4/5a repair by rescoping a shared Hopcroft-Tarjan routine
  (`dynamic_repair::find_bicc_scoped`) to just the touched SCC, not a textbook incremental biconnectivity
  algorithm — a deliberate trade-off since the higher layers already bound the affected region.
- New `dynamic_repair.h`: iterative (non-recursive) scoped Tarjan SCC and Hopcroft-Tarjan BiCC/AP,
  shared by both initial construction and incremental repair so there's one algorithm, not two.
- Unified `insertEdge`/`deleteEdge` as the two public entry points (replacing the old
  `updateEdgeDirected`/`updateEdgeUndirected` pair that could desynchronize the directed/undirected
  views); a dynamic edge overlay (`extra_*_adj` / `removed_*`) lets every traversal see inserted/deleted
  edges without needing to resize the immutable base CSR.

**Bugs found and fixed along the way** (surfaced by a new differential test, see below — none of these
were caught before because nothing checked incremental results against a from-scratch rebuild):
- `scc_core.cpp`: `nodeToScc[i]` was set to the *raw* `scc_id[]` label (a trim counter, the elephant-SCC
  sentinel, or a mice-pivot vertex id) instead of a compact index into `sccs`, so `sccs[nodeToScc[v]]`
  didn't correspond to `v`'s actual component. Combined with the label-collision bug from the code
  review, this caused genuinely wrong SCC output, not just a slower recomputation.
- `wcc_detection()` (`includes/wcc_cpu/wcc_core.cpp`) and `scc_detection()`
  (`includes/scc_cpu/scc_core.cpp`) both **miscount components even on the 18-vertex example_graph**
  (1 WCC instead of 2; 8 SCCs instead of 7 — vertices 2 and 9 have edges in both directions and should
  merge two 3-cycles into one 6-vertex SCC, but didn't). Root-causing those parallel pipelines is a
  separate, larger effort; `GIR::buildWCC()`/`buildSCC_DAG()` now use simple, unambiguously-correct
  sequential algorithms (plain union-find; `dynamic_repair::find_scc_scoped` over the whole graph) for
  construction instead. This trades away some construction-time parallelism for correctness — flagging
  as a known follow-up if profiling shows it matters on the large datasets.
- `txt_to_bin_undirected.cpp`'s never-symmetrizes-edges bug (from the code review) means the on-disk
  `undirected_*.bin` files aren't reliably symmetric. Rather than requiring the data to be regenerated,
  `GIR::buildBCT()` now symmetrizes the base undirected CSR once at load time
  (`base_undirected_adj`), so GIR's own behavior is correct regardless of that upstream file's quality.

**Testing**: `gir_test.cpp`'s `EdgeUpdateTest` (previously disabled — the only thing that would have
caught the original bugs) is rewritten as a differential test: replay random insert/delete operations,
rebuild WCC/SCC/BiCC from scratch after every single one via independent traversals, and diff against
GIR's maintained state. 3000 random operations against `data/example_graph/` pass with zero mismatches;
a separate targeted test (mirroring the paper's own worked example, plus a hand-verified intra-SCC split)
confirms Cases 1, 5b, and 5c specifically, since random testing on a small graph rarely re-creates a
genuine bridge once it's well-connected. `main.cpp` only runs this against graphs ≤10,000 vertices (it
rebuilds from scratch after every operation, so it isn't meant to run against production-scale datasets).

Known limitations, documented in `gir.h`/`gir.cpp` comments rather than hidden: Case 5c's WCC-split
repair is bounded by the whole affected WCC's size, not provably the smaller side (union-find doesn't
support cheap splits without a more advanced structure); Case-2 rank maintenance falls back to a rare
full global renumber if a local interval runs out of room, so its bound is amortized, not strict
worst-case.
