# RationalFunctionApproximation.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/dev/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8355790.svg)](https://doi.org/10.5281/zenodo.8355790)
[![codecov](https://codecov.io/gh/complexvariables/RationalFunctionApproximation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/complexvariables/RationalFunctionApproximation.jl)

[![Paper](https://proceedings.juliacon.org/papers/10.21105/jcon.00208/status.svg)](https://doi.org/10.21105/jcon.00208)

![logo](logo-sm.png)

This Julia package adaptively computes rational approximations (i.e., ratios of polynomials) for functions on intervals and other domains in the complex plane.

The [documentation](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/) includes a walkthrough showing off the main capabilities of the package.

## Quadratic Thiele approximation

This fork adds quadratic Thiele continued fractions (QTCF) to RFA v0.4.1.
Select `QTCF()` explicitly, or use `symmetric` to wrap the target and domain:

```julia
using RationalFunctionApproximation, ComplexRegions

r = approximate(symmetric(abs), symmetric(Segment(-1, 1));
    stagnation=200, max_degree=200)
```

The domain wrapper pairs points about the segment's midpoint. The function
wrapper specifies the value map, defaulting to `conj`. This implementation
retains the paper experiment routines, general derivatives, and Thiele
counterparts for construction, evaluation, conversion, and arithmetic.

See the [QTCF guide](docs/src/qtcf.md) for symmetry choices and evaluation.
The [experiment runner](test/qtcf_test.jl) reproduces the paper's
two interval comparisons and eight Schwarz-function examples:

```sh
julia --project=test -e 'using Pkg; Pkg.instantiate()'
julia --project=test test/qtcf_test.jl
julia --project=test test/qtcf_discrete_test.jl
```

Run with Julia 1.11 or later from this directory. Figures and numerical data
are written to `test/qtcf_output/figures/` and `test/qtcf_output/numbers/`.
The discrete runner writes the same figures and CSV formats to
`test/qtcf_discrete_output/`. Each case uses one fixed shared grid: 5,000
uniform parameter points, all continuum support nodes, and logarithmic
parameter neighbors. The real and imaginary intervals use exactly rotated
copies of the same grid. Both runners also write an `index.html` figure gallery.
Presets and curve definitions are in [qtcf_parameters.jl](test/qtcf_parameters.jl);
all figure, CSV, and HTML generation is in [qtcf_output.jl](test/qtcf_output.jl).

## Related work

* The [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/) package provides rational functions, but not in a way related to function approximation.
* The [BaryRational](https://juliahub.com/ui/Packages/General/BaryRational) package implements the original (fully discrete) version of the AAA algorithm, as well as Floater–Hormann rational interpolation.
* The [ApproxFun](https://juliaapproximation.github.io/ApproxFun.jl/stable) package provides 1D and multidimensional function approximation using Chebyshev polynomials and Fourier series. It also has extensive functionality for manipulating the approximations and for solving differential equations.
* There is an [ApproxFunRational](https://github.com/tomtrogdon/ApproxFunRational.jl) package, but it is undocumented.
