# RationalFunctionApproximation.jl

Julia routines for rational approximation with AAA, Thiele continued fractions,
and quadratic Thiele continued fractions (QTCF).

```julia
using RationalFunctionApproximation, ComplexRegions

r = approximate(symmetric(abs), symmetric(Segment(-1, 1));
    stagnation=200, max_degree=200)
```

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
Generated results are not included in the repository.

[qtcf_parameters.jl](test/qtcf_parameters.jl) contains the interval and Schwarz
presets, each shared across continuum and discrete modes. All fits use Float64.
The discrete grids contain 5,000 uniform parameter samples, continuum support
nodes, logarithmic neighbors, and reflected partners. The real and imaginary
interval grids are exact rotations of one another.

[qtcf_test.jl](test/qtcf_test.jl) defines the curves and fits;
[qtcf_output.jl](test/qtcf_output.jl) produces all figures and numerical tables.
