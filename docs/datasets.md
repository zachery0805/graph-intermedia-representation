# Datasets

All benchmark graphs live under `data/` (kept on disk, **not** tracked in git except the toy example —
see the repo `.gitignore`). Most are stored as CSR binary (four files: `fw_begin.bin` / `fw_adjacent.bin`
for the forward graph, `bw_begin.bin` / `bw_adjacent.bin` for the transpose, plus an `undirected_*`
triple for the symmetrized version).

| Dataset | Path | Size on disk | Notes |
|---|---|---|---|
| example_graph | `data/example_graph/` | 68K | 12-vertex toy graph, tracked in git as a test fixture |
| web-Google | `data/web-Google_sym/` | 967M | Web crawl graph, ~875K vertices / ~5M edges |
| soc-LiveJournal | `data/livejournal/` | 2.3G | Has original `soc-LiveJournal1.txt` edge list; canonical LiveJournal copy (see below) |
| RoadUSA | `data/RoadUSA_sym/` | 5.0G | US road network |
| grid_1000_100000 | `data/grid_1000_100000/` | 14G | Synthetic grid graph; also bundles its own txt<->bin conversion tools |
| twitter | `data/twitter/` | 104G | Large-scale Twitter follow graph (Kwak et al. 2010) — the same source DZiG's paper cites as its "Twitter (TW)" dataset (~1.47B edges / ~41.7M vertices) |

**Total:** ~126GB. This is why `data/` is git-ignored wholesale (only `example_graph/` is whitelisted) —
regenerate or re-copy the rest from wherever they were originally sourced.

## 2026-07-27 cleanup

Deleted to reclaim disk space (server-wide disk was at 100%, 16GB free — see CHANGELOG.md):
`data/soc-LiveJournal/` (1.3G, a second/older LiveJournal conversion, duplicate of `data/livejournal/`
which is kept as canonical), `data/CHEM_5_sym.bin` (433M, chemistry graph, unrelated to the SCC/BCC
paper line of work), `data/GeoLifeNoScale_2.bin` (1.1G, GPS-trajectory graph, likewise unrelated), and
`data/com-youtube_sym.bin` (89M, unused by the current benchmark suite). Also deleted
`data/twitter/twitter.txt` (25G raw edge-list text — the `.bin` CSR files converted from it are kept and
are what `GIR`/`bench_main.cpp` actually load at runtime; the `.txt` would need to be re-sourced if the
`.bin` files ever need to be regenerated from scratch).

Freed ~28GB total (16GB → 43GB available). This is a small dent, not a fix: the disk being full is
mostly other users'/other projects' usage on a shared server (this project's `data/` was only ~153GB out
of 6.6TB used), so adding a Friendster-scale dataset (~150-250GB estimated) is on hold pending a real
storage solution (quota increase or external storage), not something on-disk cleanup here can solve.
