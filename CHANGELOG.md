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
