# Experiments

Results and analysis, kept separate from source code so the two evolve independently and so this
directory can be pointed to directly when writing up the paper.

- `results/` — raw run output (logs, CSVs). `results/baselines/` holds prior runs of the `standard/iCC/`
  baseline suite (SCC/WCC), carried over from before the reorganization.
- `analysis/` — scripts/notebooks that turn `results/` into tables and plots. Empty for now.

See [`../docs/experiment-plan.md`](../docs/experiment-plan.md) for what's planned and
[`../docs/datasets.md`](../docs/datasets.md) for the datasets these runs are against.
