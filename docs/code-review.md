# Code Review — 2026-07-21

> **Update (2026-07-22):** the dynamic edge-update API and the SCC-DAG population bug described below
> have been fixed — see `CHANGELOG.md`'s "Rewrote GIR's dynamic edge-update (SPU) to match the paper's
> design" entry for what changed and how it was validated. The `txt_to_bin_undirected.cpp` /
> `scc_core.cpp` label-collision bugs are still present in those files as originally described (GIR now
> works around their effects internally rather than requiring them to be fixed first — see the same
> changelog entry) and the rest of this document (KSP/Yen/PeeK family, graph_converter's other findings,
> etc.) is still an open list, not yet acted on.

Deep correctness audit of `gir_constructor/` (core system) and `graph_converter/` (conversion tools),
done before the codebase is used to produce paper results. Five review passes (core GIR class, KSP/Yen/PeeK
family, Delta-Stepping/SSSP + I/O, scc_cpu/wcc_cpu baselines, graph_converter) plus a targeted adversarial
verification pass on the highest-impact claims (marked **[VERIFIED]** below — re-read directly against the
source, not just taken from the review pass).

**Bottom line:** the static graph representations and I/O are mostly solid, but **the dynamic
edge-update API (`updateEdgeDirected`/`updateEdgeUndirected` and everything under it) is essentially
non-functional**, **the SCC-DAG is never actually populated** (`dag.bin` always records 0 neighbors per
SCC, `getDAG()` always returns empty), and there's a **confirmed data-corruption bug in the graph
converters** (`txt_to_bin_undirected.cpp` never symmetrizes edges, `txt_to_bin_fw.cpp` writes to the wrong
files entirely). None of this was caught before because the only test exercising the update path is
commented out in `main.cpp`, and `wcc_cpu`'s standalone binary doesn't even compile against its own header.

Severity legend: 🔴 corrupts results silently · 🟡 crashes / fails loudly or is clearly incomplete ·
⚪ dead code / cosmetic, no output impact.

---

## 1. Core GIR class & graph representations (`gir_constructor/{gir,graph,graph_undirected,bicc_bfs,block_cut_tree,dag_graph,main}.{h,cpp}`)

### 🔴 `updateWCCDelete` is a no-op disguised as a traversal — **[VERIFIED]**
`gir.cpp:290-318`. `visited` is populated by unconditionally inserting *every* node already in
`wccs[wcc_id]` — the `std::stack` pushed with `u` is never popped/processed. So `new_component` is always
empty: deleting a bridge edge that should split a WCC into two components never does.
**Fix:** do a real BFS/DFS from `u` over the remaining edges (excluding the deleted one); anything not
reached becomes the new component.

### 🔴 `updateSCCDelete` can never fire — **[VERIFIED]**
`gir.cpp:321-336`. Guard is `if (sccs[scc_id].size() <= 1)`, but by construction `nodeToScc[u]==nodeToScc[v]`
at the call site means that SCC has ≥2 members — the guard is never true. Deleting an edge that was the
sole cycle-closing edge of an SCC should potentially split it; nothing here does that.
**Fix:** recompute strong connectivity of the SCC's induced subgraph after edge removal (local
Tarjan/Kosaraju re-run) and split `sccs[scc_id]` accordingly.

### 🔴 `isNewCycleFormed` traverses the wrong structure — **[VERIFIED bug exists; independently re-read]**
`gir.cpp:339-356`. Traverses `sccs[nodeToScc[curr]]` (the *current SCC bucket*) instead of real graph
adjacency (`g_directed->fw_beg_pos`/`fw_csr`). Called specifically when `u`/`v` are in *different* SCCs —
exactly the case this traversal can never bridge, since it only ever walks within `v`'s existing bucket.
**Fix:** BFS/DFS over actual directed edges from `v`, looking for `u`.

### 🔴 `updateSCCInsert` over-merges on new-vertex insertion
`gir.cpp:414-445`. Copies WCC-insert semantics (correct for *weak* connectivity) onto SCC membership: a
single directed edge `u→v` where `v` is new gets `v` merged into `u`'s SCC — but one directed edge does not
imply mutual reachability. A pure sink `v` should form its own singleton SCC.
**Fix:** create a new singleton SCC for an untracked endpoint; only merge across a confirmed cycle.

### 🔴 `updateBCCInsertInternal` never merges blocks
`gir.cpp:448-468`. Only checks "is `u`/`v` already in *some* BCC"; if not, creates a trivial `{node}`
singleton. Never merges components when the new edge creates a cycle spanning two existing blocks, never
touches `ap[]`/`BlockCutTree`. For endpoints already in some block, this is a complete no-op.
**Fix:** implement real BCC-insertion maintenance (cycle detection across existing blocks/articulation
points, merge on cycle).

### 🔴 `updateBCCDeleteInternal` can cascade-corrupt unrelated blocks
`gir.cpp:359-379`. Comment claims "delete nodes with degree ≤ 1" but there's no degree check — `curr` is
unconditionally erased from *every* `bicc_components` entry containing it, and "neighbors to re-examine"
are pulled from co-membership in those BCCs (not real graph adjacency), which can cascade arbitrarily far
through shared articulation points.
**Fix:** only touch the block(s) actually invalidated by the deleted edge; derive re-examination
candidates from real adjacency in `g_undirected`.

### 🔴 `updateBCCInsertCross`/`updateBCCDeleteCross` — declared, never implemented
`gir.h:83-84`. No definitions anywhere; only the `*Internal` variants are ever called. This is exactly the
harder "edge crosses between blocks" case — currently entirely unsupported.

### 🔴 SCC-DAG is never actually built — **[VERIFIED]**
`gir.cpp:142-143,197-266,474`. `buildSCC_DAG()` calls `generate_and_save_dag(...)` (the only DAG-writing
path actually invoked); inside it, the loop that inserts into `dag_adj_set[u_super]` is entirely commented
out (lines ~221-241), so `dag_adj_set[i]` is empty for every `i` and `dag.bin` records `neighbor_count = 0`
for every SCC, regardless of the real graph. Separately, `buildDagAdjList()` (the *only* function that
populates the `dag_adj_list` member returned by `getDAG()`) is never called — its call site is commented
out with the note "暂时注释掉，因为 gir.h 中的声明需要两个参数" (temporarily commented out, the
declaration needs two params). Confirmed via direct read: `getDAG()` returns `dag_adj_list`, which is
default-constructed empty and populated nowhere else in the codebase. **Any paper claim about the SCC-DAG
structure or `getDAG()` output is currently not backed by real computation.**
**Fix:** fix `buildDagAdjList`'s signature mismatch and call it from `buildSCC_DAG()` (it already has
working, correctly-locked parallel logic — just wire it in), and/or restore the commented-out population
loop in `generate_and_save_dag` (iterate vertex IDs directly rather than the `unordered_map`, to keep it
OpenMP-parallelizable).

### 🔴 `dag_graph.h` is a fully dead, unimplemented class
No `.cpp`, not in the Makefile, not `#include`d anywhere. A parallel, more complete design for
SCC-DAG+BCC+WCC that was abandoned mid-way. Would fail to link if ever instantiated.
**Fix:** implement and wire in, or delete to avoid misleading future readers.

### 🟡 Isolated (degree-0) vertices vanish from the Block-Cut Tree — **[VERIFIED]**
`bicc_bfs.h:301-306` (confirmed by direct read). For an isolated vertex, `num_bicc`/`size_map` are updated
but nothing is pushed to `bicc_components` — every other code path does push. `BlockCutTree::build` only
builds from `bicc_components`, so isolated vertices get no BCT representation at all.
**Fix:** push `{root}` onto `bicc_components` in this branch too.

### 🟡 Unchecked `std::map::operator[]` reads risk silent corruption / OOB
`gir.cpp` — `nodeToWcc[u]`/`nodeToScc[u]` used as existence checks in `updateWCCDelete`/`updateSCCDelete`/
`isNewCycleFormed`. `operator[]` auto-inserts `0` for missing keys; if a vertex never went through an
insert-update on this instance, this silently fabricates `nodeToScc[u]=0` and can index `sccs[0]` out of
bounds if `sccs` is empty/small.
**Fix:** use `.count()`/`.find()` for existence checks; handle "not tracked" explicitly.

### ⚪ Assorted lower-severity / dead-code items (details in original review, not repeated in full here)
- `graph.cpp`: fw temp buffers (`tmp_beg_pos`/`tmp_csr`) leaked before being reassigned for bw section —
  add `delete[]` before reuse.
- `graph.cpp`/`graph_undirected.cpp`: vestigial `index_tt`→`index_t` copy is currently a no-op (both are
  `int`) — dead conversion path, safe to simplify.
- `graph.h`: `gen_src()`, `groupby()` empty no-ops; `graph_load()` duplicates constructor logic and is
  unused — dead API surface, delete or implement.
- `bicc_bfs.h`: `get_size_result` only ever writes the first histogram bucket instead of iterating all of
  `size_map` — output under-reports BCC size distribution.
- `bicc_bfs.h`: local `#define DEBUG 0`/`OUTPUT_TIME 0` after `util.h` already defined them to `1` means
  `make DEBUG=1` can never actually enable this header's debug blocks.
- `bicc_bfs.h`: unused globals/dead code (`time_index`, `check_root`, and the `checked`/`low`/`par` arrays
  threaded through every function but never read).
- `gir.cpp`: hardcoded, inconsistent thread counts (6 vs 16) and unexplained `alpha/beta/gamma/theta`
  tuning constants — hurts reproducibility, should come from config/CLI.
- `gir.cpp`: `buildBCT()` leaks `g_undirected`/`ap`/`bc_tree` if called twice (no `delete` before
  reassignment).
- `main.cpp`: `#include "gir_test.cpp"` (including a `.cpp` directly) — fragile, ODR risk if ever added to
  the Makefile's `SRCS`.
- `main.cpp:69-72`: the only test of the update API (`EdgeUpdateTest::testUpdateEdgeDirected`) is
  commented out — directly explains why the bugs above went unnoticed.

---

## 2. KSP / Yen / PeeK family (`gir_constructor/includes/{KSP,Yen,OptYen,PeeKAdaptive*,PeeKEdgeSwap,PeeKRegeneration,PeeKStatusArray,KSPGraph,MaxHeapNoDuplicate}.h`)

**Context: this whole family is currently orphaned** — nothing in the compiled build (`main.cpp` →
`gir.cpp`/`block_cut_tree.h`) references `Yen`/`OptYen`/any `PeeK*` class; the only reference anywhere is
`legacy/bc_tree.h`, which isn't compiled either. That's *why* these bugs were never caught — but they'd all
resurface immediately if this family is revived for the paper.

### 🔴 No guard for "destination unreachable" → unsigned underflow → crash
`Yen.h:105,168`; `OptYen.h:179,235,806,862`. The first candidate path is seeded and popped unconditionally.
If source/destination are disconnected, `SsspTree::get_shortest_path()` returns an empty vector, and the
deviation loop's bound `paths[kp].p.size() - 1` (unsigned) underflows to a huge value — first-iteration
out-of-bounds access on an empty vector.
**Fix:** after seeding, check `if (paths[0].length == INFINITE_DISTANCE || paths[0].p.empty()) return paths;`
before the deviation loop, at all four call sites.

### 🔴 Garbage `Path` candidates can be inserted and later crash the same way
`Yen.h:388-464`; `OptYen.h:422-535,1027-1140`. Parallel candidate buffers default-construct empty/`INFINITE_
DISTANCE` `Path` slots; unfilled slots (common when true candidate count < k) get bulk-inserted with no
validity filter, and can later be selected via `get_best()`.
**Fix:** filter `candidates_buffer` before the bulk insert, or add the validity check inside
`Candidates::insert(Path&&)` itself.

### 🟡 `PeeKAdaptiveWithEdgeSwap`'s "adaptive" branch is dead — always takes the edge-swap path
`PeeKAdaptiveWithEdgeSwap.h:81-93`. The `if (k_bound_nodes_num > adaptive_ratio * ...)` switch is commented
out, so the entire "regenerate compact graph" branch (fully written, lines 95-150) is unreachable. The
sibling `PeeKAdaptiveWithEdgeSwapParallel.h` has the *same* `if` intact — clear evidence this was
disabled during debugging and never restored.
**Fix:** uncomment the guard to match the Parallel sibling.

### 🟡 Leaked `CSRGraph` on the graph-regeneration path
`PeeKRegeneration.h:290`, `PeeKAdaptiveWithStatusArray.h:305`, `PeeKAdaptiveWithEdgeSwap.h:306`,
`PeeKAdaptiveWithEdgeSwapParallel.h:440`. `->destroy()` frees internal arrays but not the `CSRGraph` object
itself — leaked once per query on that code path.
**Fix:** `delete compact_g.get_csr_graph();` after `destroy()`, or own it via `unique_ptr`.

### ⚪ Data race on a plain `bool` (`is_loop`) in the parallel edge-swap variant
`PeeKAdaptiveWithEdgeSwapParallel.h:238-244` — `volatile bool` read/written across threads with no
atomics; UB in principle, and the "early stop" it implies doesn't actually stop other threads anyway.
**Fix:** `std::atomic<bool>`, or drop it (path lengths here are small, serial scan is fine).

### ⚪ Copy-paste artifacts
- `PeeKAdaptiveWithEdgeSwap.h` and `PeeKAdaptiveWithEdgeSwapParallel.h` define the **same class name**
  (`PeeKAdaptiveWithEdgeSwap`) — hard redefinition error if both are ever included together. Rename the
  Parallel file's class to match its filename.
- `PeeKStatusArray.h:90` — the only variant with an unconditional (uncommented) debug
  `optyen.print_paths();`; every sibling has it commented out.
- `MaxHeapNoDuplicate.h`'s two headline classes are entirely unused dead code; the free function `cmp()`
  is a non-inline header-defined function — ODR landmine if ever pulled into two translation units.
- All 6 `PeeK*.h` files' `normalSSSP()` permanently removes edges from `orig` via `remove_node`/
  `edge_is_removed` and never calls `restore_graph()` — harmless today (each caller constructs a fresh
  object per query) but a landmine if `compute()` is ever called twice on one instance.

---

## 3. Delta-Stepping SSSP / graph I/O (`gir_constructor/includes/{DeltaStepping*,BucketListStatic,RequestListStatic,Request,Partition,SsspTree,BasicGraph,GraphRW,NodeList,Node,Profiling,Statistics}.h`, `gir_test.cpp`)

### 🔴 Data race on `std::vector<bool> node_is_removed`
`DeltaStepping.h:305`, `DeltaSteppingStatic.h:126-131,191-196`. `vector<bool>` is bit-packed; concurrent
writes to different logical bits that share a storage word race even when each *node* is only written by
its owning thread — two nodes owned by different threads can land in the same word at partition
boundaries. `#pragma omp flush` orders visibility, not atomicity of the underlying word RMW.
**Fix:** use `std::vector<char>` or `std::atomic<bool>`/`vector<std::atomic<bool>>` instead.

### 🔴 `delta_option` parameter is completely ignored; delta silently becomes 0 for unweighted graphs
`GraphRW.h:33,149-151,181,235-237` + `BasicGraph.h:210-213`. Neither `read_graph` nor `read_graph_bct`
reference their own `delta_option` parameter — delta is unconditionally computed as
`avg_weight/avg_degree`, and `GraphRW::get_delta()` (which implements the documented `-D`/`-BF`/explicit
selection) is never called anywhere. Worse: for `weighted=false` instantiations, `total_weight` stays 0,
so `delta` computes to exactly `0`, and `compute_bucket_size()` then divides `heaviest_weight/delta` —
`inf`, then UB on the `double→unsigned int` conversion.
**Fix:** actually call `GraphRW::get_delta(delta_option, heaviest_weight)`; guard the unweighted case with
a fixed sane delta instead of falling through the averaging formula with `total_weight==0`.

### 🔴 `CSRGraph*` leaked from every `read_graph`/`read_graph_bct`/`get_reverse_graph` call
`BasicGraph.h` (no destructor at all) + `GraphRW.h:75,189,257`. A comment on `get_reverse_graph` claims
"the object gets deleted in the destructor of the original graph" — no such destructor exists.
**Fix:** give `BasicGraph` a real destructor that deletes `csr`; fix copy semantics accordingly.

### 🟡 Dead code path calls `SsspTree` methods that don't exist
`DeltaSteppingStatic.h:53,69` call `sssp_tree.define_parent_path(...)`/`update_branch()` — not declared
anywhere in `SsspTree.h`, confirmed via repo-wide grep. This overload is never actually instantiated by any
current call site, so it doesn't break the build today, but would fail to compile the moment anyone uses
it.

### 🟡 `partition_edges_by_weight`'s unweighted fast path is permanently dead
`GraphRW.h:345-397`. Every call site invokes it without explicit template args, so `weighted` always
resolves to its `true` default — the `!weighted` branch can never run regardless of the graph's actual
weighted-ness.
**Fix:** forward the real `weighted` argument from `read_graph<directed,weighted>` at each call site.

### 🟡 `PartitionDynamic::get_thread_id` is a stub returning 0
`Partition.h:115-126`, author's own comment admits "does this even work?". Not currently wired into
`DeltaStepping<...>`'s template parameters anywhere — dormant, but would cause every node to route to
thread 0 (assert failures or a real data race) if ever plugged in.

### ⚪ Statistics / test coverage gaps
- `DeltaSteppingStatic.h:78` (PROF builds): `touched` initialized to `1` instead of `0` — systematic +1
  bias in every `nodes_touched_in_ds` measurement.
- `Profiling.h`: ~16 collected fields, only 3 surfaced via `return_stats()`; `print_times()` fully
  commented out; `nodes_touched_in_ds` has an explicit `// TODO does not work in pd mode`.
- `Partition.h:78-82`: `rand_nodes` is never actually randomized despite the name/todo — contiguous
  chunking only.
- `NodeList.h:13-18`: `clear_buffer()` is a no-op in release builds (only asserts an invariant, doesn't
  enforce it).
- **`gir_test.cpp` has zero coverage of the entire Delta-Stepping/SSSP/GraphRW subsystem** — grepped for
  every relevant symbol, zero hits. `main.cpp`'s only SSSP test reference is commented out and calls a
  function (`ParallelSSSP::testParallelSSSP`) that doesn't even exist in `gir_test.cpp`.
- `gir_test.cpp:331-334`: fabricated "per-phase timing" — four printed lines are the *same* single
  `bfs_time` measurement multiplied/offset by unexplained magic constants (`*50`, `*20-1.654`, `*10`,
  `*20-15`), labeled as if they were four independently-measured phases. Also: `testParallelBFS` never
  validates BFS correctness against a reference, only counts/prints `visited_count`.

---

## 4. SCC/WCC baselines (`gir_constructor/includes/{scc_cpu,wcc_cpu}/`)

**Context:** unlike KSP/PeeK, this subsystem is load-bearing — `gir.cpp`'s `buildWCC()`/`buildSCC_DAG()`
call into it, and the README documents `batch_test*.sh` scripts for it.

### 🔴 SCC label ID-space collision can silently merge unrelated SCCs — **[VERIFIED]**
`scc_cpu/scc_core.cpp:38` + `fw_bw.h`. Confirmed by direct read: `atomic_scc_label` starts at `2` and
increments for every trimmed vertex/group; the largest ("elephant") SCC is hardcoded to `scc_id = 1`; and
`mice_fw_bw` labels small residual SCCs using **the pivot vertex's own numeric ID**
(`scc_id[v] = v`/`scc_id[w] = v`). All three draw from the same `[0, vert_count)` integer space with no
reserved offset — a trim-assigned ID can coincide with a real vertex ID that later becomes a mice pivot,
silently merging two unrelated SCCs in `scc_core.cpp`'s final grouping. Vertex `0` is also a broken
sentinel: `scc_id[vert]==0` means "unresolved" *and* is a legitimate pivot label if vertex 0 is ever chosen.
**Fix:** reserve disjoint ID ranges (e.g. offset the atomic counter by `vert_count+2`, offset mice-pivot
labels by `2*vert_count`), or tag IDs with a "kind" discriminant instead of relying on raw integer equality.

### 🟡 `wcc_main.cpp` does not compile against the current `wcc_core.cpp`/`wcc_common.h` — **[VERIFIED]**
Confirmed by direct read: `wcc_common.h` declares `void wcc_detection(graph*, index_t, vertex_t*)` (3
args, `extern "C"`); `wcc_main.cpp:43-47` calls it with 5 arguments of unrelated types
(`avg_time, result_file_wcc, result_file_time` where a single `vertex_t*` is expected), and never even
allocates a `wcc_id` array. **The standalone `wcc_cpu` executable described in the README cannot currently
build at all** — this looks like `wcc_core.cpp`'s API was refactored for the `gir.cpp` integration without
updating `wcc_main.cpp`. Any prior claim that `wcc_cpu`'s batch tests ran against the current code should
be treated as unverified.
**Fix:** update `wcc_main.cpp` to allocate a `wcc_id` array and call the 3-argument signature; confirm it
actually builds before trusting any of its batch-test output.

### 🟡 Dead loop guard defeats the `theta` convergence parameter
`scc_cpu/scc_core.cpp:419-457`. An unconditional `break;` at the end of the trim-repeat loop means it
always runs exactly one iteration regardless of the `theta`-based convergence check — compare with the
structurally identical (and correct) loop at lines 284-326. Doesn't produce wrong final SCCs
(`mice_fw_bw` still resolves everything) but silently makes `theta` partly inert here.

### 🟡 Unbounded thread-local frontier buffer, no bounds check
`scc_cpu/fw_bw.h:130` (and the wcc analog). `thread_queue` is sized by a heuristic
(`vert_count/thread_count * 5` or `*7`) with **no bounds check** on writes — a high-degree hub vertex or
benign frontier-duplication race can overrun this heap buffer.
**Fix:** bound writes explicitly, or just size to `vert_count` (safe upper bound).

### ⚪ Broken/misleading statistics and dead code
- `load_print.cpp`'s `get_scc_result`: classifies by sentinel values (`-1/-2/-3/-4`) that the trim
  functions no longer produce (refactored to the shared atomic counter) — the printed size-1/2/3
  breakdown is **always 0**; everything falls into the "others" bucket. Any paper number quoting this
  breakdown is currently wrong.
- `scc_core.cpp`: uninitialized `time_gfq` read into `avg_time[14]` (UB; currently masked because the
  corresponding printf is commented out).
- `scc_core.cpp`: local `#define DEBUG 0` comes *after* headers that already expanded `DEBUG` to `1` at
  the point of use in `scc_common.h` — a stray debug line prints on every run despite the apparent intent
  to suppress it.
- `wcc_cpu/{trim_1_gfq,trim_2_3,fw_bw,color_propagation}.h`: fully dead, `#include`s commented out in
  `wcc_core.cpp` — ~2500 lines of SCC-flavored code copy-pasted in and never used for WCC.
  `wcc_core.cpp` itself also carries 5 alternate WCC algorithms as dead code (`wcc_hybrid` is the only one
  actually called).
- `scc_cpu/color_propagation.h`: ~350 lines of an alternate, never-invoked SCC-identification pipeline
  (only the final `coloring_wcc` function in the file is used).
- `scc_cpu/trim_2_3.h`: `_from_graph` variants (whole-graph scan) are dead; only `_from_fq` variants run.
- `scc_main.cpp:21-23`: `alpha`/`beta`/`gamma` parsed with `atof` then truncated to `const int` — silently
  drops fractional CLI values.
- SUSPECTED (not confirmed with a concrete failing input): hardcoded path-compression depth caps
  (`c_depth < 100`/`< 50`) in the coloring phases could leave compression incomplete on very deep union
  chains; unsynchronized concurrent reads/writes of `color[]`/`wcc_id[]` during label propagation are a
  real data race by the memory model, but look like an intentional lock-free convergent design (both
  writers converge to the same value) rather than an output-corrupting bug — flagging for awareness, not
  as confirmed-wrong.

---

## 5. `graph_converter/` — the tools that produce every binary `gir_constructor` loads

This is the highest-stakes area: bugs here silently corrupt the input to every experiment.

### 🔴 `txt_to_bin_undirected.cpp` never symmetrizes edges — **[VERIFIED]**
Confirmed by direct read (`edge_list_to_csr_bin/txt_to_bin_undirected.cpp:143-188`): the reverse-edge
insertion (`adj[begin[v_id]+degree[v_id]] = index; degree[index_a]++;`) is commented out at **both** the
counting pass and the fill pass. Only the forward direction (`index → v_id`) is ever stored. This tool
writes to the exact filenames (`undirected_begin.bin`/`undirected_adjacent.bin`) that `gir_constructor`
loads for `g_undirected` — **every "undirected" graph produced by this tool is actually still directed**,
silently. Any WCC/BCC computation that assumes true symmetric adjacency is wrong for any vertex whose
edges are all incoming in the original directed listing.
**Fix:** uncomment the reverse-edge lines (both the counting-pass `degree[index_a]++` and the fill-pass
`adj[begin[v_id]+degree[v_id]]`/`head[...]` writes) — this exact correct logic already exists, working, in
`bin_to_bin/bin_to_bin_undirected.cpp` and can be copied from there.

### 🔴 `txt_to_bin_fw.cpp` writes to the wrong output files entirely — **[VERIFIED]**
Confirmed via `diff`: `txt_to_bin_fw.cpp` is byte-for-byte identical to `txt_to_bin_undirected.cpp`. It
opens `undirected_adjacent.bin`/`undirected_head.bin`/`undirected_out_degree.bin`/`undirected_begin.bin`
— the `fw_*.bin` filenames are present only as dead commented-out lines. **Running the tool named/built
specifically to produce the forward CSR (`fw_begin.bin`/`fw_adjacent.bin`) never creates those files** —
`gir_constructor` will either fail to open them (hard `exit(-1)`, the safe case) or silently load a stale
copy from a previous run (the dangerous case).
**Fix:** restore this file to build the forward CSR under the correct names (model it on
`txt_to_bin_fw_int.cpp`, which is correct, just with wider integer types) — delete the
`undirected_*`-writing version or rename/remove the file to prevent confusion.

### 🔴 8-byte vs 4-byte integer width mismatch between converters and `gir_constructor`'s reader
`edge_list_to_csr_bin/txt_to_bin_bw.cpp` and `tuple_text_to_bin.cpp` use 8-byte
`vertex_t`/`index_t` (`unsigned long long`/`long int`); `gir_constructor/graph.cpp` hardcodes 4-byte
`int` reads (`index_tt`/`vertex_tt`) and computes `vert_count` from file size ÷ element size. If
`bw_begin.bin` ever comes from `txt_to_bin_bw` (8-byte) instead of its `_int` sibling (correct, 4-byte),
`vert_count` computes to roughly double the truth and every offset/adjacency value gets read at the wrong
byte boundary — either an `assert` abort (safe) or, if sizes coincidentally line up, silently scrambled
backward-adjacency data.
**Fix:** delete the 8-byte variants or rename them clearly so they can't be mistaken for
`gir_constructor`-compatible output; standardize on the `_int` (4-byte) tools.

### 🟡 `graph_clean.cpp`'s dedup only works on input pre-sorted by source vertex
`edge_list_to_csr_bin/graph_clean.cpp:22-59`. The `set` used for duplicate detection is cleared every time
the source vertex changes — non-adjacent duplicate edges (e.g. `1 2` ... `2 3` ... `1 2` again) are never
recognized as duplicates. No sortedness precondition is checked or documented anywhere.
**Fix:** require+validate sorted input, or replace the streaming dedup with a full
`unordered_set<pair<int,int>>` / sort-then-unique pass.

### 🟡 No CSR-builder performs its own dedup/self-loop stripping
None of `txt_to_bin_{fw,bw}_int.cpp`/`txt_to_bin_bw.cpp`/`txt_to_bin_undirected.cpp` check for duplicate
pairs or self-loops — they rely entirely on the separate (and, per above, buggy-for-unsorted-input)
`graph_clean` step. Running a CSR builder directly on raw data (an easy mistake, nothing enforces the
two-step pipeline) lets self-loops/duplicates flow straight into `fw_adjacent.bin`/`bw_adjacent.bin`.

### 🟡 `bin_parser.cpp` is an explicit, self-documented stub
`bin_to_bin/bin_parser.cpp:16-48`. Comments literally say "临时硬编码，仅用于演示文件结构" (temporarily
hardcoded, just to demonstrate the file structure) — hardcodes `V_count`/`E_count` instead of reading them
from the file; the real `fs.read()` calls are commented out. Produces meaningless output for any graph
other than the one specific shape it was hacked together for.

### 🔴 `g500_bin_to_edge_list.py` reads int64 vertex IDs as int32 — every edge corrupted
Confirmed by cross-referencing `batch_g500.py`'s invocation of the standard Graph500
`generator_test_mpi` (which emits 8-byte `packed_edge{long v0, v1}` records — same struct shape appears,
commented out, in this codebase's own `tuple_text_to_bin.cpp:17-20`) against `g500_bin_to_edge_list.py`'s
`chunk_size = 4` / `struct.unpack('i', ...)`. Splitting each real 8-byte ID into two 4-byte reads means
every output "edge" has its second endpoint forced to 0 for realistic (< 2^32) vertex counts — **any
Graph500 dataset run through this script produces a star graph centered on vertex 0**, not the real graph.
**Fix:** `chunk_size = 8`, `struct.unpack('q', ...)`, matching what `rmat_bin_to_edge_list.py` already
does correctly for its own (also 8-byte) input format.

### ⚪ Assorted lower-severity items
- `graph_converter/graph.h`: dead, non-compilable duplicate of `gir_constructor/graph.h` — stray clutter.
- `tuple_text_to_bin.cpp`: the *only* tool with correct symmetrization logic, but writes to
  non-`gir_constructor`-matching filenames and uses the mismatched 8-byte width — looks like the correct
  logic ended up in the wrong file during some copy-paste history.
- `bin_to_bin_bw_int.cpp`/`bin_to_bin_fw_int.cpp`: despite living under `bin_to_bin/` and being named like
  conversion tools, they perform no transformation — they verify+byte-copy the same data under new
  filenames.
- `bin_to_bin_undirected.cpp`: otherwise correct symmetrization, but doesn't dedup if the input directed
  graph already contains both `u→v` and `v→u` as separate edges (common for "mutual friendship"-style
  datasets) — would double-count those edges. Impact depends on which datasets route through this vs. the
  (broken) `txt_to_bin_undirected`.
- `convert_graphs.sh`: hardcoded absolute path tied to one user's home directory.
- `edge_list_to_matrix_market.py`: never writes the mandatory `%%MatrixMarket` header line — output isn't
  a valid MM file, standard readers (`scipy.io.mmread`) will reject it.
- **Python 2/3 split**: `directed_edge_list_to_undirected{,_back}.py` (`from sets import Set`),
  `batch_g500.py` (`.iteritems()`), and several scripts' bare `print x, y` statements are all Python-2-only
  syntax — will not run at all under Python 3. `graph_remove_size_zero.py` is already Python-3-clean,
  so the codebase is inconsistent and needs a single porting pass.
- `directed_edge_list_to_undirected.py` vs `_back.py`: two files doing the same job with different
  internal approaches, no indication which is canonical; the non-`_back` one keeps self-loops
  (inconsistent with `graph_clean.cpp`'s convention).
- `unsorted_edgelist_to_sorted.py`: hardcoded single-space delimiter, while the C++ tools tolerate
  space/tab/newline — will break on tab-separated input.
- `graph_correct_edge.py`: no filename filtering before parsing (any non-matching file aborts the whole
  batch), unquoted shell `mv` command.
- `split_vertex_group.py`: brittle fixed-position filename parsing, no length validation.
- `edge_list_to_adj.py`: expects a bare `v_count e_count` header line while
  `directed_edge_list_to_undirected*.py` write a `%`-prefixed header — chaining these two scripts breaks.

---

## Priority order for fixing (before relying on results for the paper)

1. **`txt_to_bin_undirected.cpp` / `txt_to_bin_fw.cpp`** (graph_converter) — these corrupt the *input
   data itself*; nothing downstream can be trusted until these are fixed and existing converted datasets
   under `data/` are re-generated.
2. **SCC label collision** (`scc_cpu`) — silently merges unrelated SCCs; affects any reported SCC/WCC
   count or the (currently also broken) DAG.
3. **DAG never populated** (`gir.cpp`) — fix before reporting anything about the SCC-condensed DAG.
4. **The dynamic update path** (`updateWCCDelete`/`updateSCCDelete`/`isNewCycleFormed`/
   `updateSCCInsert`/`updateBCCInsertInternal`/`updateBCCDeleteInternal`) — if the paper's contribution
   involves dynamic/incremental maintenance, this is currently the weakest part of the whole codebase and
   needs the most work; write real tests for it first (re-enable/fix `main.cpp`'s commented-out
   `EdgeUpdateTest` calls) so future fixes can be checked against something.
5. **`wcc_main.cpp` build break** — confirm what `wcc_cpu`'s batch-test results in
   `experiments/results/baselines/wcc/` were actually run against, since the standalone binary doesn't
   currently compile.
6. Everything else (KSP/PeeK family, dead code, statistics/profiling gaps) — lower urgency since it's
   either orphaned (not in the current build) or cosmetic, but worth cleaning up before it's revived for
   later experiments.

Findings not marked **[VERIFIED]** above come from a single thorough review pass with concrete file/line
citations but were not independently re-read by a second pass — treat them as high-confidence leads, not
already-double-checked facts, and re-verify against the current code before citing them in the paper.
