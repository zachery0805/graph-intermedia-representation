# Datasets

All benchmark graphs live under `data/` (kept on disk, **not** tracked in git except the toy example —
see the repo `.gitignore`). Most are stored as CSR binary (four files: `fw_begin.bin` / `fw_adjacent.bin`
for the forward graph, `bw_begin.bin` / `bw_adjacent.bin` for the transpose, plus an `undirected_*`
triple for the symmetrized version).

| Dataset | Path | Size on disk | Notes |
|---|---|---|---|
| example_graph | `data/example_graph/` | 68K | 12-vertex toy graph, tracked in git as a test fixture |
| com-youtube | `data/com-youtube_sym.bin/` | 89M | YouTube social network, ~1.1M vertices / ~3M edges |
| web-Google | `data/web-Google_sym/` | 967M | Web crawl graph, ~875K vertices / ~5M edges |
| CHEM_5 | `data/CHEM_5_sym.bin/` | 433M | Chemistry graph |
| GeoLifeNoScale_2 | `data/GeoLifeNoScale_2.bin/` | 1.1G | GeoLife trajectory-derived graph |
| soc-LiveJournal (v1) | `data/livejournal/` | 2.3G | Has original `soc-LiveJournal1.txt` edge list |
| soc-LiveJournal (v2) | `data/soc-LiveJournal/` | 1.3G | ~4.8M vertices / ~69M edges — check which of the two LiveJournal copies is canonical before running experiments, they look like separate conversions of the same source |
| RoadUSA | `data/RoadUSA_sym/` | 5.0G | US road network |
| grid_1000_100000 | `data/grid_1000_100000/` | 14G | Synthetic grid graph; also bundles its own txt<->bin conversion tools |
| twitter | `data/twitter/` | 129G | Large-scale Twitter follow graph |

**Total:** ~153GB. This is why `data/` is git-ignored wholesale (only `example_graph/` is whitelisted) —
regenerate or re-copy the rest from wherever they were originally sourced.

## Open item

`data/livejournal/` and `data/soc-LiveJournal/` both look like LiveJournal conversions but aren't
identical (different file sets, no `.txt` in the second one). Worth resolving which one is authoritative
before it causes confusion in an experiment writeup — flagging here rather than silently picking one.
