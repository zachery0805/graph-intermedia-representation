# Legacy

Pre-`gir_constructor` prototypes, superseded by the code in `gir_constructor/` and kept here only for
historical reference (they predate the point where `gir.cpp`, `graph.cpp`, etc. were consolidated into
`gir_constructor/`). Do not build on these — extend `gir_constructor/` instead.

- `bcc.cpp`, `bc_tree.h`, `bicc_bfs.h` — early Block-Cut-Tree / biconnected-components prototype.
- `findSCC.cpp` — early standalone SCC prototype.
- `gir_construction.cpp`, `gir.cpp`, `graph.cpp`, `graph.h` — early `GIR` prototype, since replaced by
  `gir_constructor/gir.cpp` and `gir_constructor/graph.cpp` (different implementation, not just a copy).
- `util.h`, `wtime.h` — shared utility headers used by the above.
- `bash.sh` — build script for the early prototype (`g++ gir.cpp -o gir`).
- `graph.py`, `scc.py` — early Python prototyping scripts.

## `flagged/`

Files quarantined during the 2026-07 repository reorganization for the repo owner to review and delete
once confirmed unneeded — not reviewed/deleted automatically since they're unrelated to graph algorithms:

- `check_permissions.sh` — a generic sysadmin script (iterates `/etc/passwd`, checks per-user
  permissions/sudo access) found under `standard/iCC/tempol/`. Looks like it was copied in by mistake;
  has nothing to do with GIR or the iCC baseline suite.
