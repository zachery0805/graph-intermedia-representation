# Changelog

Milestone-level summary of major changes. Day-to-day detail lives in `git log`; this file is for
skimming "what happened when" while writing up the paper.

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
