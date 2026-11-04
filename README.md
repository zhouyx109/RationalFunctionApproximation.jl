# RationalFunctionApproximation.jl

Julia routines for rational approximation with AAA, Thiele continued fractions,
and quadratic Thiele continued fractions (QTCF).

```julia
using RationalFunctionApproximation, ComplexRegions

r = approximate(symmetric(abs), symmetric(Segment(-1, 1));
    tol=1e-11, stagnation=50, max_degree=200,
    initial_refinement=15, refinement=4)
```

Symmetric QTCF starts at the midpoint by default. Supply `fixed_points=[z₁, z₂, …]`
to use each fixed point as a separate scalar startup update, in the supplied
order. Subsequent updates use symmetric pairs and skip other fixed points.
For a sampled domain, the midpoint or supplied startup points must be present
in the grid. `initial_node` can override the default midpoint.

## Paper experiments

Run from this directory with Julia 1.11 or later:

```sh
julia --project=test -e 'using Pkg; Pkg.instantiate()'
julia --project=test test/qtcf_test.jl
julia --project=test test/qtcf_discrete_test.jl
```

The runners reproduce the two interval comparisons and eight Schwarz-data
examples. Each writes four figures in PDF and PNG, numerical CSV files, and an
`index.html` gallery to `test/qtcf_output/` or `test/qtcf_discrete_output/`.
Figures use the paper's fixed pole windows and residue color scale.
Generated results are not included in the repository.

[qtcf_parameters.jl](test/qtcf_parameters.jl) contains the interval and Schwarz
presets, each shared across continuum and discrete modes. All fits use Float64.
Discrete interval fits preserve the original TCF experiment's ordered union of
positive and negative logarithmic samples `2.^(-52:1//10:-1)`, 10,001 uniform
points, and logarithmic samples shifted toward the left endpoint. The 11,509
points are neither sorted nor completed by reflection. The imaginary grid is
exactly `im` times the real grid. Both methods use the same supplied points.
Interval validation uses the fixed grid specified in the paper in both fitting modes.

Discrete Schwarz grids contain 5,000 uniform parameter subdivisions, the
interpolation nodes from every recorded continuum stage, and logarithmic
neighbors at the offsets in `DISCRETE_LOG_EXPONENTS`. Parameters are wrapped on closed curves
and clamped on open curves. Reflected partners complete each grid exactly,
and all methods use the same grid for each case.
In both experiment modes QTCF starts with one scalar block at the parameter midpoint,
then use symmetric pairs. The five-lobed parameter interval is
`[-pi/2, 3pi/2]`, so its midpoint is fixed by the prescribed reflection.
Independent validation retains fixed points.

QTCF stores polynomial numerator and denominator blocks. Evaluation, derivatives,
poles and residues use these same stored coefficients; singleton storage delegates
to Thiele's kernels. Values and derivatives use the backward block recurrence.
Selected-node scalar and array values at common zeros use local polynomial factor cancellation.
Pole extraction defaults to the unreduced denominator pencil;
`poles(r; cancel_common=true)` removes confirmed common factors using the same
pencil eigenvalues. Residues use the shared Thiele routine.
Degree queries report nominal bounds propagated from the block sizes, as TCF
does for one-point blocks; actual degrees may be lower.
The experiments use `OneDiv` evaluation, `OneDiv` new-weight computation for TCF,
and `Classic` new-weight computation for QTCF. The parameter file sets
the tolerance, degree limit, refinement, and stagnation limit for each family.
At the current degree limit of 200, the update limits are 400 nodes for TCF,
200 blocks for QTCF, and 200 updates for AAA.
QTCF supports generic floating-point data; the experiments use Float64.
Separate continuum and discrete loops use an adaptive curve mesh or
a supplied fixed grid. Both select the largest remaining residual,
require finite target values, final reduced values and block coefficients,
and save every completed block. Error checks include the remaining samples
and previously selected nodes. The sampled target maximum is updated at every
block. Continuum refinement decreases once per block and resamples all
subintervals. Both modes return the best acceptable recorded approximant.
Continuum TCF starts with 15 interior samples per subinterval.
Both use OneDiv array sweeps with scalar retries for NaN values. QTCF also
recovers nonfinite values at selected common zeros. Scalar evaluation retains
the shared underflow protection. Independent validation evaluates each point
with the scalar routine.

Numerical tables include validation and sampled errors, stopping reasons,
returned records, extraction status, and applied fitting and sampling settings.
`run_metadata.csv` records the source commit, working-tree status, source and
dependency-file hashes, package versions, and runtime precision. Pole tables
retain every extracted pole and residue, including nonfinite results.
Independent interval validation uses the paper's fixed set `X1 ∪ X2 ∪ (-X2) ∪
(X2-1)`, where `X1` has 10,001 uniform points on `[-1,1]` and
`X2 = 2.^(-0.1 .* (10:1000))`. The imaginary interval uses `im .* X`.
Schwarz validation uses 1,800 uniform parameter points on every curve,
counting the common endpoint once for closed curves. These fixed test sets
are shared by all methods, all recorded interval stages, and both fitting
modes; they do not depend on the selected interpolation nodes.
Convergence curves end at the smallest validation error and show denominator
degrees from 0 to 100. All computed history records remain in the CSVs.

[qtcf_test.jl](test/qtcf_test.jl) defines the curves and fits;
[qtcf_output.jl](test/qtcf_output.jl) produces all figures and numerical tables.
