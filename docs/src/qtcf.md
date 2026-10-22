```@meta
CurrentModule = RationalFunctionApproximation
```

# Quadratic Thiele continued fractions

[`QuadraticThiele`](@ref), also named `QTCF`, builds rational approximations
with one- or two-point continued-fraction blocks. It supports the continuum
and discrete paper experiments, general derivatives, and the usual Thiele
construction and evaluation interfaces:

```julia
using RationalFunctionApproximation, ComplexRegions

r = approximate(symmetric(abs), symmetric(Segment(-1, 1)))
z = collect(range(-1, 1, 101))
r = approximate(symmetric(abs), symmetric(z; reflection=x -> -x))
```

Results support the usual `get_function`, [`status`](@ref),
[`isconverged`](@ref), [`get_history`](@ref), [`rewind`](@ref),
`domain`, `test_points`, and [`check`](@ref) operations.

## Point and value symmetry

[`symmetric`](@ref) wraps either a domain or a function. A domain wrapper
specifies a point map; a function wrapper specifies a value map. Together
they declare

```math
f(\rho(z)) = \sigma(f(z)).
```

```julia
r = approximate(
    symmetric(abs),
    symmetric(Segment(-1, 1));
    stagnation=200, max_degree=200,
)
```

For a `Segment`, the default point map reflects about its midpoint:
`z -> 2c - z`. The default function-value map is `conj`. In this example,
`abs(-x) == conj(abs(x))`. Either wrapper selects QTCF by default.

Supply `reflection` to choose another map:

```julia
r = approximate(symmetric(cos; reflection=identity), symmetric(Segment(-1, 1)))
r = approximate(symmetric(sin; reflection=y -> -y), symmetric(Segment(-1, 1)))
```

Other domain types require an explicit point reflection. [`Symmetric`](@ref)
is the underlying wrapper type. A plain target with a wrapped domain evaluates
its partner values directly; wrapping the target supplies the value map.
The target's value relation is
checked at sampled reflected points. These checks do not prove symmetry
throughout the domain.

A discrete grid must contain every reflected partner. The discrete fit uses
only the supplied node-value data; it does not generate off-grid partners.

Continuum QTCF locates each reflected point on the curve and inserts both selected
node parameters into the mesh before refinement. It uses analytic inverses or
verified parameter correspondences where available, otherwise the geometry
library's numerical locator. Failure to locate a partner is reported rather than
assigning it the source parameter. The initial partition is `[0, 1]` for a curve
or `[0, length(domain)]` for a multi-curve path. The first block selects
`point(domain, 0)` and its reflection only if distinct; its selected parameters
then refine the partition. Other fixed points and curve breakpoints are not
preloaded. Test samples lie inside the mesh subintervals, as in TCF; unselected
endpoints are excluded, while a supplied symmetric startup pair is retained.
QTCF has no golden-section routine.
The point map is assumed to be an involution preserving the approximation domain.
The stored `result.path` contains this partition; QTCF also completes reflected
samples during iteration. Generic `check` samples the partition, while the paper
experiments use independent validation grids.

## Construction and evaluation

Direct construction defaults to singleton blocks. Use `block_sizes` for pairs:

```julia
z = [-1.0, 1.0, 0.0]
g = QTCF(z, abs.(z); block_sizes=[2, 1])
RFA.add_node!(g, [-0.5, 0.5], [0.5, 0.5])
```

`weights(g)` returns one constant or affine coefficient vector per block.
`nodes`, `values`, `degrees`, `degree`, `copy`, `isreal`, numeric conversion,
and scalar arithmetic follow the Thiele interface. `QTCF{Float32}(z, y)` and
`convert(BigFloat, g)` select coefficient precision. Supplied weights can be
passed as the third constructor argument, with the same `block_sizes`.
Without a symmetry map, `approximate(f, domain, QTCF())` starts with one node,
then adds the two largest-residual candidates per block.

Obtain the interpolant from the approximation result:

```julia
g = get_function(r)
g(0.2)
evaluate(g, 0.2)
evaluate(g, [0.2, 0.4], Classic())
RFA.evaluate!(zeros(2), g, [0.2, 0.4])
derivative(g, 4)(0.2)
derivative(g, [0, 1, 2, 4])(0.2)
poles(g)
roots(g)
h = convert(PartialFractions, g)
locations, strengths = residues(g)
```

Scalar evaluation uses the one-division recurrence with TCF's underflow safeguard
and backward-evaluation fallback. Singleton representations delegate to TCF.
[`derivative`](@ref) returns an evaluator for one order or a vector of orders.
For mixed blocks, recovery of all requested derivatives reuses one denominator
reciprocal; singleton representations use Thiele's derivative routines.
Numerator roots and denominator roots use the same block matrix-pencil routine; residue computation uses the
numerator and denominator derivative ratio, with contour integration as a fallback.
Partial-fraction conversion retains the polynomial part and requires simple poles.
Exact common factors at selected nodes are cancelled during pole extraction.
Selected-node evaluation also resolves removable zero-over-zero values.

Construction reuses RFA's result and history types. Continuum fits reuse sampling
and reduced-value buffers; residue extraction reuses its derivative workspace.

## Approximation controls

The usual `tol`, `float_type`, `allowed`, and `stagnation` keywords apply.
Continuum QTCF also accepts `refinement`, `initial_refinement`, and
`max_degree`; `max_iter` defaults to `max_degree` and counts blocks, not
denominator degree. For sampled data, use `max_iter`.

`allowed=true` permits all poles; a predicate tests individual poles.
For continuum domains, `allowed=:strict` checks distance from a curve or
exclusion from a region.
As in TCF, check `isconverged(r)` before treating an unsuccessful return as
meeting the tolerance and pole constraints.

New blocks require finite target values, final reduced values, and stored
coefficients. Every convergent is saved. Error checks include previously selected
nodes, while only unselected points can enter a new block.

## Paper experiments

From the package directory, with Julia 1.11 or later:

```sh
julia --project=test -e 'using Pkg; Pkg.instantiate()'
julia --project=test test/qtcf_test.jl
julia --project=test test/qtcf_discrete_test.jl
```

The runner wraps both functions and domains with `symmetric` for every QTCF
case. It runs only the paper's two interval comparisons and eight
Schwarz-function examples: 28 fits, four figures (PDF and PNG), and eleven CSV
files. Rerunning replaces the files in `test/qtcf_output/figures/` and
`test/qtcf_output/numbers/`. The paper runner is separate from the upstream
package tests in `test/runtests.jl`.

`qtcf_test.jl` runs the fits and builds the shared discrete grids.
`qtcf_parameters.jl` stores presets, curve definitions, and display settings;
`qtcf_output.jl` generates all figures, CSV files, and the HTML gallery.
The discrete entry point reuses the same three files.
Including the runner loads its definitions; `run_experiments()` starts the fits.

The discrete runner first constructs the continuum approximants, then fits
the same data on one fixed shared grid per case. Following the stable-evaluation
sampling pattern, the grid combines parameters `j/5000` for `j=0:4999`, every
method's continuum support nodes, and parameter offsets `±10^k` for
`k=-16:0.1:-3` around those nodes. Neighbors are mapped through the curve;
closed curves wrap at the endpoint and open curves retain points in `[0,1]`.
The grid includes exact reflected partners for every method. Exact duplicates
are removed without merging distinct logarithmic neighbors. The two interval
cases share the union of their continuum support and neighbor samples, with
the imaginary grid obtained by multiplying the real grid by `im`.

Discrete fits use the discrete presets in `qtcf_parameters.jl`, without adaptive
refinement. Error checks use the whole fixed grid, including selected nodes,
and the requested tolerance relative to the target scale. For equal-value symmetry
after a fixed point, discrete selection also uses the preceding paired error.

The plotting and validation code is shared with the continuum
runner. Outputs go to `test/qtcf_discrete_output/`, with
`numbers/discrete_grid_metadata.csv` recording grid sizes and stopping reasons.
Each output directory contains an `index.html` figure gallery.
