# Experiment Plan

Source: originally `gir_constructor/tasks/experiments.txt`, expanded here. This is the running list of
experiments intended for the GIR paper. Update this file as experiments are added, run, or dropped —
raw output goes in `experiments/results/`, analysis/plots go in `experiments/analysis/`.

## Planned experiments

1. **Throughput** — sustained update/query throughput of GIR under continuous edge insertions/deletions.
2. **Concurrent BFS** — parallel BFS performance on the directed/undirected CSR graph representations.
3. **Concurrent SSSP** — Delta-Stepping SSSP performance (see `gir_constructor/includes/DeltaStepping*.h`).
4. **Concurrent BC** — betweenness centrality (or comparable connectivity-based) benchmark.
5. **Batch update** — applying edge updates in batches vs. one at a time; measure maintenance cost of
   SCC/WCC/Block-Cut-Tree decompositions.
6. **Single update** — per-edge update latency (insert/delete), same decompositions as above.

## Baselines for comparison

- `standard/iCC/` — in-house baseline suite (SCC/WCC/BiCC/bridge detection, CPU and GPU variants).
  Prior run outputs live in `experiments/results/baselines/`.
- `~/projects/ppopp16` — "Multi-Core On-The-Fly SCC Decomposition" (external reference implementation).
- `~/projects/Parallel_SSSP_Algorithm` — reference parallel dynamic SSSP update algorithm.

## Datasets

See [`datasets.md`](datasets.md) for the full list and where each lives on disk.

## Status

Not yet tracked per-experiment — add a table here (dataset x experiment x status) once runs start,
so it's easy to see what's done when writing the paper.
