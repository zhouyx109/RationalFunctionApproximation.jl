# Continuum QTCF experiments

The main comparison runs 27 Float64 fits across 11 targets. It uses midpoint-start QTCF and endpoint-start Paired TCF on every target, and native square TCF on the five point-symmetry targets. Paired TCF greedily chooses a node and then forces its reflected partner as a scalar insertion. Native square TCF runs plain TCF on the folded target in `u = z^2`. All methods use Classic construction and OneDiv evaluation.

## Rerun

Use a Git checkout containing commit `e19e600b4d0bca9ed12d6e6fa6369923ee6329b6`, Julia 1.11.5, Python 3.9 or later, and the usual `bash`, `patch`, and `tar` commands. If using a shallow clone, fetch the pinned commit. Install the Python dependencies in your active environment and run from the repository root.

```sh
python3 -m pip install pypdf==6.14.2 sympy==1.14.0
./experiments/qtcf_continuum/run.sh ./qtcf-results
```

Supply a new directory, or omit it to create a fresh results directory in the current directory. The launcher applies `rfa_experiment.patch` to a temporary package snapshot and instantiates its pinned Julia test environment. Set `QTCF_INSTANTIATE=0` if those dependencies are already installed. It uses one Julia and one BLAS thread. Wall-clock times depend on the machine. If the run stops, the launcher prints the locations of the logs and work files.

## What the files do

- `run.jl` fits the targets and records errors, timing, selected nodes and weights. It times three fits after warming up and runs the diagnostic trace separately.
- `check_matched_midpoint.jl` repeats the six Group L Paired TCF fits at the QTCF midpoint. Its table is separate from the 27-fit comparison.
- `claim_1.py` replays scalar TCF in exact arithmetic at the recorded QTCF nodes. For the three odd targets, insertion breaks down at the first pair. The checks at two and four old pairs use separately seeded scalar prefixes. At one and three old pairs the recorded reduced type cannot be produced by finite scalar insertion on those nodes.
- `stress_3_3.jl` reruns QTCF on the two absolute-value intervals at refinement 3/3. `replay_interval_reductions.py` compares the 28/9 interval reduction traces with exact arithmetic.
- The degree and pole scripts independently check the stored rational functions. The figure and table scripts use those checks to build the comparison.

## Outputs

- `summary.csv` reports errors and fit times in seconds. The errors use original coordinates. Divide `error` by `validation_target / 4e-11` to obtain the plotted relative error.
- `histories.csv` and `reduced_history_degrees.csv` record the histories by physical node count and exact reduced denominator degree. `diagnostics.csv` and `reductions.csv` give attempted weights and outcomes.
- `poles.csv` records raw pencil results. `verified_curve_poles.csv` records the high-precision pole and residue checks used for the curve plots. The `interval_certification/` folder holds the exact interval checks.
- `matched_midpoint.csv` compares Group L errors, nodes, reduced degrees and stored denominator symmetry at the same startup.
- `seeded_claim.csv`, `same_order_replay.csv`, and `stress_3_3*.csv` record the breakdown, same-order, and stress checks.
- `qtcf_required_comparison_figures.pdf` and `group_*_rows.tex` contain the regenerated plots and table rows. The individual logs and `provenance.txt` record the run settings and audit steps.
