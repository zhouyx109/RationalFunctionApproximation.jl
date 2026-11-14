# -*- coding: utf-8 -*-
# %% [markdown]
# # Quadratic Thiele continued fractions
#
# ## 1. Types and accessors
#
# `_QTCFBlocks` stores blocks; `QuadraticThiele` uses Thiele storage for singletons.
# `nodes`, `values`, `degrees`, and `degree` support RFA results and display.
# `convert` changes numeric type; `copy` saves snapshots; `isreal` checks coefficients.
# `getproperty` builds continuants for polynomial queries.

# %%
import Base: values

export QuadraticThiele, QTCF, Symmetric, symmetric

# %% _QTCFBlocks — block storage
# QTCF-only storage: denominator coefficients are vectors, one per block.
struct _QTCFBlocks{S}
    pairs::Vector{Tuple{S,S}}
    bpolys::Vector{Vector{S}}
    apolys::Vector{Vector{S}}
    nodes::Vector{S}
    vals::Vector{S}
end

# Curve and fixed-grid data used to form symmetry blocks.
struct _QTCFCurveSampler{D,F,R}
    domain::D
    target::F
    parameter_reflection::R
end

struct _QTCFDiscreteSampler{Y,Z}
    values::Y
    points::Z
    remaining::Vector{Int}
end

# %% QuadraticThiele — approximant type
"""
    QuadraticThiele()
    QuadraticThiele(nodes, values, [weights]; block_sizes, method=default_weight_method())

Quadratic Thiele representation, aliased as `QTCF`. Pass the empty selector to
`approximate`, or construct from node-value data. `block_sizes` contains ones
and twos (default: all ones); optional weights are constant or affine coefficient
vectors, one per block. Singleton storage delegates to Thiele.
Continuants are built lazily from the stored blocks.
"""
mutable struct QuadraticThiele{T,S} <: AbstractRationalInterpolant{T,S}
    data::Union{Thiele{T,S},_QTCFBlocks{S}}
    polynomials::Union{Nothing,Tuple{Vector{S},Vector{S}}}

    function QuadraticThiele{T,S}(data, polynomials) where {T<:AbstractFloat,S<:RealComplex{T}}
        return new{T,S}(data, polynomials)
    end
end

# Construct from scalar or block storage.
QuadraticThiele(data::Thiele{T,S}, polynomials) where {T,S} =
    QuadraticThiele{T,S}(data, polynomials)
QuadraticThiele(data::_QTCFBlocks{S}, polynomials) where {S} =
    QuadraticThiele{real_type(S),S}(data, polynomials)

# %% Base.convert
"""
    convert(F, r::QuadraticThiele)

Convert the numeric type of a QTCF.
"""
function Base.convert(::Type{F}, r::QuadraticThiele{T,S}) where {F<:AbstractFloat,T,S}
    d = getfield(r, :data)
    d isa Thiele && return QuadraticThiele(convert(F, d), nothing)
    V = S <: Real ? F : complex(F)
    data = _QTCFBlocks{V}(Tuple{V,V}.(d.pairs), [V.(b) for b in d.bpolys],
        [V.(a) for a in d.apolys], V.(d.nodes), V.(d.vals))
    return QuadraticThiele{F,V}(data, nothing)
end

# %% nodes
"Return the selected QTCF interpolation nodes in stored block order."
nodes(r::QuadraticThiele) = r.nodes

# %% Base.values
"Return the target values at the QTCF interpolation nodes."
values(r::QuadraticThiele) = r.vals

# %% weights
"Return the expanded constant or affine denominator-block coefficients in block order."
weights(r::QuadraticThiele) = r.bpolys

# %% degrees
"Return the nominal QTCF numerator and denominator degree bounds."
function degrees(r::QuadraticThiele)
    d = getfield(r, :data)
    d isa Thiele && return degrees(d)

    pold, qold = 0, typemin(Int)
    p, q = (first(d.pairs)[1] == first(d.pairs)[2] ? 0 : 1), 0
    for j in 2:length(d.pairs)
        dj = d.pairs[j][1] == d.pairs[j][2] ? 1 : 2
        previous = d.pairs[j-1][1] == d.pairs[j-1][2] ? 1 : 2
        p, q, pold, qold = max(dj - 1 + p, previous + pold),
            max(dj - 1 + q, previous + qold), p, q
    end
    return p, q
end

# %% degree
"Return the nominal QTCF denominator degree bound."
degree(r::QuadraticThiele) = last(degrees(r))

# %% Base.isreal
# Check the stored nodes and weights, as for Thiele.
function Base.isreal(r::QuadraticThiele)
    data = getfield(r, :data)
    return data isa Thiele ? isreal(data) : isreal(data.nodes) && all(isreal, data.bpolys)
end

# %% Base.copy
# Copy the stored data, blocks, and any cached continuants.
function Base.copy(r::QuadraticThiele{T,S}) where {T,S}
    d = getfield(r, :data)
    data = d isa Thiele ? copy(d) : _QTCFBlocks{S}(
        copy(d.pairs), copy.(d.bpolys), copy.(d.apolys), copy(d.nodes), copy(d.vals))
    p = getfield(r, :polynomials)
    return QuadraticThiele{T,S}(data, p === nothing ? nothing : copy.(p))
end

# %% Base.getproperty
# Return stored nodes/values or build polynomial continuants.
@inline function Base.getproperty(r::QuadraticThiele{T,S}, name::Symbol) where {T,S}
    d = getfield(r, :data)
    if name === :P || name === :Q
        p = getfield(r, :polynomials)
        if p === nothing
            b, a = d isa Thiele ? (d.weights, d.nodes) : (d.bpolys, d.apolys)
            Pold, Qold = S[1], S[0]
            P, Q = isempty(b) ? (S[0], S[1]) :
                (d isa Thiele ? S[first(b)] : copy(first(b)), S[1])

            # Advance the Euler--Wallis recurrence without cancelling factors.
            for j in 2:length(b)
                bj, aj = d isa Thiele ? ((b[j],), (-a[j-1], one(S))) : (b[j], a[j-1])
                P, Q, Pold, Qold = _padd(_pmul(bj, P), _pmul(aj, Pold)),
                    _padd(_pmul(bj, Q), _pmul(aj, Qold)), P, Q
            end
            r.polynomials = p = (P, Q)
        end
        return name === :P ? first(p) : last(p)
    end
    name in (:data, :polynomials) && return getfield(r, name)
    name === :values && (name = :vals)
    name === :weights && (name = :bpolys)
    if d isa Thiele
        name === :nodes && return d.nodes
        name === :vals && return d.values
        name === :pairs && return [(z, z) for z in d.nodes]
        name === :bpolys && return [[w] for w in d.weights]
        name === :apolys && return [S[-z, 1] for z in d.nodes]
    end
    return getproperty(d, name)
end

# %% QTCF
"Alias for [`QuadraticThiele`](@ref)."
const QTCF = QuadraticThiele

# %% QuadraticThiele() — method selector
# Create an empty selector without entering a data-construction API.
QuadraticThiele() = QuadraticThiele{Float64,Float64}(Thiele(), nothing)

# %% [markdown]
# ## 2. Evaluation, polynomial support, and derivatives
#
# Scalar `evaluate` uses OneDiv/Classic; array evaluation uses OneDiv. `_evaluate_numden`
# and `_evaluate_classic` supply their block recurrences; `_evaluate_block` evaluates
# one block. `_evaluate_common_zero` cancels a selected linear factor locally.
# `_padd` and `_pmul` build continuants in `getproperty`.
# `_evaluate_block_derivs` and `_evaluate_numden_derivs` supply polynomial
# derivatives for `residues`; `_derivative!` supplies general quotient derivatives.

# %% evaluate — scalar
"""
    evaluate(r::QuadraticThiele, z, method=default_eval_method())

Evaluate one point with `OneDiv()` or backward `Classic()` evaluation.
OneDiv retries Classic at a zero denominator; scalar evaluation rescales the polynomial pair.
At a selected common zero, evaluate after canceling its simple common factor.
Infinity uses the nominal block degree bounds, as scalar Thiele uses node count.
"""
function evaluate(r::QuadraticThiele, z::Number, method::ThieleMethod=default_eval_method())
    d = getfield(r, :data)
    d isa Thiele &&
        return convert(promote_type(eltype(r), typeof(z)), evaluate(d, z, method))

    # Resolve nonfinite inputs before entering the block recurrence.
    if isinf(z)
        P, Q = r.P, r.Q
        m, n = degrees(r)
        p = m < length(P) ? P[m+1] : zero(eltype(r))
        q = n < length(Q) ? Q[n+1] : zero(eltype(r))
        return m > n ? oftype(last(P), Inf) :
            m < n ? zero(last(P)) : p / q
    elseif isnan(z)
        return oftype(z, NaN)
    end

    return _evaluate_finite_thiele(r, z, method)
end

# %% evaluate — array
# Promote stored coefficients with the query type.
function evaluate(r::QuadraticThiele, z::AbstractArray{<:Number})
    t = similar(z, promote_type(eltype(r), eltype(z)))
    return _evaluate_thiele_samples!(t, r, z, similar(t), similar(t), default_eval_method())
end

# %% evaluate! — allocate scratch
"""
    evaluate!(t, r::QuadraticThiele, z, [a, b])

Write evaluations into `t` using matching, nonoverlapping scratch arrays.
With `fallback=true`, retry nonfinite results and small intermediate values
through scalar OneDiv. Selected common zeros use scalar cancellation.
"""
function evaluate!(
    t::AbstractArray, r::QuadraticThiele, z::AbstractArray{<:Number}
    )
    return evaluate!(t, r, z, similar(t), similar(t))
end

# %% evaluate! — reuse scratch
# Fill t using caller-supplied scratch arrays a and b.
function evaluate!(
    t::AbstractArray, r::QuadraticThiele, z::AbstractArray{<:Number},
    a::AbstractArray, b::AbstractArray;
    fallback=false
    )
    d = getfield(r, :data)

    d isa Thiele && return evaluate!(t, d, z, a, b; fallback)

    retry = fallback ? fill!(similar(t, Bool), false) : nothing
    m = length(d.bpolys)
    if m == 1
        p = d.bpolys[1]
        if length(p) == 1
            t .= p[1]
        else
            t .= p[1] .+ p[2] .* z
        end
    else
        p = d.bpolys[m]
        if length(p) == 1
            a .= p[1]
        else
            a .= p[1] .+ p[2] .* z
        end
        p, pair = d.apolys[m-1], d.pairs[m-1]
        if length(p) == 2
            b .= p[2] .* (z .- pair[1])
        else
            b .= (p[3] .* (z .- pair[1])) .* (z .- pair[2])
        end
        _mark_thiele_underflow!(retry, a, b)
        @inbounds for k in m-1:-1:2
            p = d.bpolys[k]
            if length(p) == 1
                b .= p[1] .* a .+ b
            else
                b .= (p[1] .+ p[2] .* z) .* a .+ b
            end
            p, pair = d.apolys[k-1], d.pairs[k-1]
            if length(p) == 2
                a .*= p[2] .* (z .- pair[1])
            else
                a .*= (p[3] .* (z .- pair[1])) .* (z .- pair[2])
            end
            a, b = b, a
            _mark_thiele_underflow!(retry, a, b)
        end
        p = d.bpolys[1]
        if length(p) == 1
            t .= p[1] .+ b ./ a
        else
            t .= (p[1] .+ p[2] .* z) .+ b ./ a
        end
    end

    return _finish_thiele_array!(t, r, z, retry)
end

# %% _evaluate_classic
# Evaluate by the backward continued-fraction recurrence.
function _evaluate_classic(r::QuadraticThiele, z::Number)
    d = getfield(r, :data)
    d isa Thiele && return _evaluate_classic(d, z)
    m = length(d.bpolys)
    u = _evaluate_block(d.bpolys[m], d.pairs[m], z)
    for k in m-1:-1:1
        u = _evaluate_block(d.bpolys[k], d.pairs[k], z) +
            _evaluate_block(d.apolys[k], d.pairs[k], z; numerator=true) / u
    end
    return u
end

# %% _evaluate_numden
# Return the backward numerator/denominator pair with common scaling.
function _evaluate_numden(r::QuadraticThiele, z::Number; underflow_guard=true)
    @assert isfinite(z)
    d, S = getfield(r, :data), promote_type(eltype(r), typeof(z))
    if d isa Thiele
        return _evaluate_numden(d, z; underflow_guard)
    end

    # Seed the final block and carry the pair backward without division.
    b, a, pairs = d.bpolys, d.apolys, d.pairs
    m = length(b)
    A = S(_evaluate_block(b[m], pairs[m], z))
    m == 1 && return A, one(S)
    B = S(_evaluate_block(a[m-1], pairs[m-1], z; numerator=true))
    @inbounds for k in m-1:-1:1
        t = _evaluate_block(b[k], pairs[k], z) * A + B
        B = k == 1 ? A : A * _evaluate_block(a[k-1], pairs[k-1], z; numerator=true)
        A = t
        if k > 1 && underflow_guard && guard(A)
            A *= 1e100
            B *= 1e100
        end
    end
    return A, B
end

# %% _evaluate_common_zero
# Divide a selected common factor out of the local lower-tail polynomials.
function _evaluate_common_zero(r::QuadraticThiele, z)
    pairs, b, a = r.pairs, r.bpolys, r.apolys
    j = findfirst(pair -> z == first(pair) || z == last(pair), pairs)
    m = length(b)
    (j === nothing || j == m) && return nothing
    S = promote_type(eltype(r), typeof(z))
    bjet = k -> _evaluate_block_derivs(b[k], pairs[k], z)
    ajet = k -> _evaluate_block_derivs(a[k], pairs[k], z; numerator=true)
    return _cancel_thiele_factor(S, j, m, bjet, ajet)
end

# %% _evaluate_block
# Evaluate a constant/affine denominator or a factored numerator block at z.
@inline function _evaluate_block(p, pair, z; numerator=false)
    if numerator
        return length(p) == 2 ? p[2] * (z - pair[1]) :
            p[3] * (z - pair[1]) * (z - pair[2])
    end
    return length(p) == 1 ? only(p) : p[1] + p[2] * z
end

# %% _padd
# Return the trimmed sum of two polynomial coefficient vectors.
function _padd(p, q)
    n = max(length(p), length(q))
    r = zeros(promote_type(eltype(p), eltype(q)), n)
    @views r[1:length(p)] .+= p
    @views r[1:length(q)] .+= q
    return resize!(r, something(findlast(!iszero, r), min(1, length(r))))
end

# %% _pmul
# Return the trimmed product of two polynomial coefficient vectors.
function _pmul(p, q)
    r = zeros(promote_type(eltype(p), eltype(q)), length(p) + length(q) - 1)
    @inbounds for i in eachindex(p), j in eachindex(q)
        r[i+j-1] += p[i] * q[j]
    end
    return resize!(r, something(findlast(!iszero, r), min(1, length(r))))
end

# %% _evaluate_block_derivs
# Return a block value and its first two derivatives.
@inline function _evaluate_block_derivs(p, pair, z; numerator=false)
    S = promote_type(eltype(p), typeof(z))
    if !numerator
        return S(_evaluate_block(p, pair, z)), length(p) > 1 ? S(p[2]) : zero(S), zero(S)
    elseif length(p) == 2
        c, z1 = S(p[2]), S(pair[1])
        return c * (z - z1), c, zero(S)
    end
    c, z1, z2 = S(p[3]), S(pair[1]), S(pair[2])
    return c * (z - z1) * (z - z2), c * ((z - z1) + (z - z2)), 2c
end

# %% _thiele_jets! — block derivatives
# Fill numerator and denominator derivatives for the shared quotient routine.
function _thiele_jets!(z, r::QuadraticThiele, order, A, B)
    n = order + 1
    d, S = getfield(r, :data), promote_type(eltype(r), typeof(z))
    d isa Thiele && return _thiele_jets!(z, d, order, A, B)
    m = length(d.bpolys)
    fill!(view(A, 1:n), 0)
    fill!(view(B, 1:n), 0)
    bjet = _evaluate_block_derivs(d.bpolys[m], d.pairs[m], z)
    ajet = m > 1 ? _evaluate_block_derivs(d.apolys[m-1], d.pairs[m-1], z;
        numerator=true) : (one(S), zero(S), zero(S))
    for q in 0:min(order, 2)
        A[q+1], B[q+1] = bjet[q+1], ajet[q+1]
    end
    for k in m-1:-1:1
        block = d.bpolys[k]
        bjet = _evaluate_block_derivs(block, d.pairs[k], z)
        a = k > 1 ? d.apolys[k-1] : nothing
        ajet = k > 1 ? _evaluate_block_derivs(a, d.pairs[k-1], z; numerator=true) : nothing
        for q in order:-1:0
            value = bjet[1] * A[q+1]
            q > 0 && length(block) > 1 && (value += q * bjet[2] * A[q])
            value += B[q+1]
            next_B = A[q+1]
            if k > 1
                next_B = ajet[1] * A[q+1]
                q > 0 && (next_B += q * ajet[2] * A[q])
                q > 1 && length(a) > 2 &&
                    (next_B += (q * (q - 1) ÷ 2) * ajet[3] * A[q-1])
            end
            A[q+1], B[q+1] = value, next_B
        end
        k > 1 && _guard_thiele_jets!(view(A, 1:n), view(B, 1:n))
    end
    return nothing
end

# %% _evaluate_numden_derivs
# Return commonly normalized numerator/denominator jets for residues.
function _evaluate_numden_derivs(r::QuadraticThiele, z::Number)
    d = getfield(r, :data)
    d isa Thiele && return _evaluate_numden_derivs(d, z)
    m = length(d.bpolys)
    a, aʹ, _ = _evaluate_block_derivs(d.bpolys[m], d.pairs[m], z)
    m == 1 && return a, one(a), aʹ, zero(a)
    b, bʹ, _ = _evaluate_block_derivs(d.apolys[m-1], d.pairs[m-1], z; numerator=true)
    a, b, aʹ, bʹ = _normalize_residue_jets(z, a, b, aʹ, bʹ)
    @inbounds for k in m-1:-1:2
        c, cʹ, _ = _evaluate_block_derivs(d.bpolys[k], d.pairs[k], z)
        e, eʹ, _ = _evaluate_block_derivs(d.apolys[k-1], d.pairs[k-1], z; numerator=true)
        t = cʹ * a + c * aʹ + bʹ
        bʹ = eʹ * a + e * aʹ
        aʹ = t
        t = c * a + b
        b = a * e
        a = t
        a, b, aʹ, bʹ = _normalize_residue_jets(z, a, b, aʹ, bʹ)
    end
    c, cʹ, _ = _evaluate_block_derivs(d.bpolys[1], d.pairs[1], z)
    return c * a + b, a, cʹ * a + c * aʹ + bʹ, aʹ
end

# %% derivative — single order
"""
    derivative(r::QuadraticThiele, order=1)
    derivative(r::QuadraticThiele, orders::AbstractVector)

Return an evaluator for the requested nonnegative derivative order(s).
Scalar storage uses Thiele's polynomial recurrence; quotient recovery uses one reciprocal.
At selected common zeros, cancel the simple common factor first. The denominator
after cancellation must be nonzero.
"""
function derivative(r::QuadraticThiele, order::Integer=1)
    f = derivative(r, [order])
    return z -> only(f(z))
end

# %% derivative — requested orders
# Return requested derivatives, reusing captured workspaces.
derivative(r::QuadraticThiele, orders::AbstractVector{<:Integer}) = _thiele_derivative(r, orders)

# %% [markdown]
# ## 3. Poles, roots, and residues
#
# `roots` and `poles` assemble the numerator and denominator matrix pencils.
# `residues` reuses RFA's residue loop with `_evaluate_numden_derivs`.

# %% poles
"""
    poles(r::QuadraticThiele; full_pencil=false, scale=false, refine=false, finite_only=true)

Return finite roots of the unreduced denominator from the manuscript pencil
`K_m = C_m + t D_m`, with one row per block after the first, solved as the generalized
eigenvalue problem `-C_m v = λ D_m v` with no scaling and no refinement.
`full_pencil` is accepted for compatibility and uses the same pencil.
`scale=true` and `refine=true` are opt-in departures from the manuscript: joint
power-of-two row and column scaling before the solve, and Newton correction after it.
Set `finite_only=false` to return every eigensolver value, including Inf and NaN.
"""
function poles(r::QuadraticThiele; full_pencil=false, scale::Bool=false, refine::Bool=false,
    finite_only=true)
    C, D = _qtcf_block_pencil(r, 2)
    isempty(C) && return eltype(r)[]
    return _thiele_pencil_values(r, -C, D; scale, refine, finite_only)
end

# %% roots
"""
    roots(r::QuadraticThiele, first_block=1; scale=false, finite_only=true)

Return finite numerator (`first_block=1`) or denominator (`first_block=2`) roots.
Both cases use the unreduced manuscript pencil with no scaling; `scale=true` opts in
to row and column scaling.
"""
function roots(r::QuadraticThiele, first_block::Integer=1; scale::Bool=false, finite_only=true)
    C, D = _qtcf_block_pencil(r, first_block)
    isempty(C) && return eltype(r)[]
    return _thiele_pencil_values(r, -C, D; scale, refine=false, finite_only)
end

# %% residues
"""
    residues(r::QuadraticThiele; check_finite=false, radius=nothing,
        full_pencil=false, scale=false, refine=false)

Return poles and residues using the shared TCF residue loop.
"""
residues(r::QuadraticThiele; check_finite=false, radius=nothing, full_pencil=false,
    scale::Bool=false, refine::Bool=false) =
    residues(r, poles(r; full_pencil, scale, refine); check_finite, radius)

# %% [markdown]
# ## 4. Data construction, node updates, and blocks
#
# Constructors and `add_node!` reuse `_new_weight` for reduced values and
# `_qtcf_append_block!` for storage.
# Both approximation methods call `add_node!` for block construction.
# `_block_line` fits denominator coefficients; `_factor_poly` forms numerators.
# `_qtcf_from_blocks` wraps the blocks returned by arithmetic operations.

# %% QuadraticThiele — typed construction
"""
    QuadraticThiele{T}(nodes, values, [weights]; kwargs...)

Construct a QTCF with floating-point type `T`.
"""
function QuadraticThiele{T}(args...; kwargs...) where {T<:AbstractFloat}
    return convert(T, QuadraticThiele(args...; kwargs...))
end

# %% QuadraticThiele — interpolation data
# Construct from node values and optional weights, using the requested block sizes.
function QuadraticThiele(
    x::AbstractVector, y::AbstractVector,
    w::Union{Nothing,AbstractVector}=nothing;
    block_sizes=nothing, method::ThieleMethod=default_weight_method()
    )
    # Check that the block layout accounts for every datum and supplied weight.
    length(x) == length(y) || throw(DimensionMismatch("node and value counts must agree"))
    sizes = block_sizes === nothing ? ones(Int, length(x)) : collect(block_sizes)
    all(d -> d isa Integer && d in (1, 2), sizes) ||
        throw(ArgumentError("block_sizes must contain only 1 or 2"))
    sum(sizes; init=0) == length(x) ||
        throw(DimensionMismatch("block_sizes must sum to the number of nodes"))
    w === nothing || length(w) == length(sizes) ||
        throw(DimensionMismatch("provide one weight per block and its block_sizes"))

    if isempty(x)
        types = (eltype(x), eltype(y), w === nothing ? eltype(y) : eltype(w))
        types = map(T -> T <: AbstractVector ? eltype(T) : T, types)
        S = promote_type(map(T -> T === Any ? Float64 : float(T), types)...)
        return QuadraticThiele{real_type(S),S}(Thiele(S[], S[], S[]), nothing)
    end

    # Normalize scalar or polynomial weights and promote all data together.
    scalar = block_sizes === nothing || all(==(1), sizes)
    w === nothing || all(v isa Number || 1 <= length(v) <= d for (v, d) in zip(w, sizes)) ||
        throw(ArgumentError("a block weight must be constant or affine for its block size"))
    b = w === nothing ? nothing : scalar ?
        [float(v isa Number ? v : only(v)) for v in w] :
        [v isa Number ? [float(v)] : float.(collect(v)) for v in w]
    xx, yy = float.(x), float.(y)
    S = promote_type(eltype(xx), eltype(yy),
        b === nothing ? eltype(xx) : scalar ? eltype(b) : promote_type(map(eltype, b)...))
    z, f = convert(Vector{S}, xx), convert(Vector{S}, yy)

    all(isfinite, z) || throw(ArgumentError("interpolation nodes must be finite"))
    length(Set(w + zero(w) for w in z)) == length(z) ||
        throw(ArgumentError("interpolation nodes must be distinct"))

    # Singleton blocks use the Thiele representation.
    if scalar
        all(isfinite, f) || throw(ArgumentError("QTCF block target values must be finite"))
        v = b === nothing ? copy(f) : convert(Vector{S}, b)
        if b === nothing
            @inbounds for i in eachindex(v)
                v[i] = _new_weight(method, view(z, 1:i-1), view(v, 1:i-1), z[i], f[i])
                isnan(v[i]) && throw(NaNException("Adding QTCF block caused NaN reduced value"))
                isfinite(v[i]) || throw(ArgumentError("QTCF block reduced values must be finite"))
            end
        end

        any(isnan, v) && throw(NaNException("Adding QTCF block caused NaN weight"))
        all(isfinite, v) || throw(ArgumentError("QTCF block weights must be finite"))
        return QuadraticThiele{real_type(S),S}(Thiele(z, f, v), nothing)
    end

    # Build mixed blocks in order, using supplied coefficients when available.
    r, k = QuadraticThiele{real_type(S),S}(Thiele(S[], S[], S[]), nothing), 1
    for (j, size) in enumerate(sizes)
        d = Int(size)
        zz, yy = view(z, k:k+d-1), view(f, k:k+d-1)
        w === nothing ? add_node!(r, zz, yy, method) : _qtcf_append_block!(r, zz, yy, b[j])
        k += d
    end
    return r
end

# %% _qtcf_from_blocks
# Wrap owned block arrays returned by arithmetic operations.
function _qtcf_from_blocks(pairs, bpolys, apolys, z::Vector{S}, y) where {S}
    data = _QTCFBlocks{S}(pairs, bpolys, apolys, z, y)
    return QuadraticThiele{real_type(S),S}(data, nothing)
end

# %% add_node! — scalar node
"""
    add_node!(r::QuadraticThiele, z, y, method=default_weight_method())

Append one node or a vector of one or two nodes and their target values.
Nodes must be finite and distinct; target and reduced values and the resulting
block weights must be finite. Previously selected nodes are not audited for attainment.
"""
function add_node!(
    r::QuadraticThiele, z_new::Number, y_new::Number,
    method::ThieleMethod=default_weight_method()
    )
    # Validate the proposed datum before reducing it.
    z, y = eltype(r)(z_new), eltype(r)(y_new)
    isfinite(z) || throw(ArgumentError("interpolation nodes must be finite"))
    isfinite(y) || throw(ArgumentError("QTCF block target values must be finite"))
    z in nodes(r) && throw(ArgumentError("interpolation nodes must be distinct"))

    h = _new_weight(method, r, z, y)
    isnan(h) && throw(NaNException("Adding QTCF block caused NaN reduced value"))
    isfinite(h) || throw(ArgumentError("QTCF block reduced values must be finite"))
    _qtcf_append_block!(r, (z,), (y,), (h,))
    return r
end

# %% add_node! — node block
# Append a vector of one or two node-value pairs and return the updated QTCF.
function add_node!(
    r::QuadraticThiele, z_new::AbstractVector,
    y_new::AbstractVector, method::ThieleMethod=default_weight_method()
    )
    # Route singleton updates through the scalar method and validate paired data.
    length(z_new) == length(y_new) ||
        throw(DimensionMismatch("node and value counts must agree"))
    length(z_new) in (1, 2) ||
        throw(ArgumentError("a QTCF update requires one or two nodes"))
    length(z_new) == 1 && return add_node!(r, only(z_new), only(y_new), method)
    z, y = eltype(r).(z_new), eltype(r).(y_new)
    all(isfinite, z) || throw(ArgumentError("interpolation nodes must be finite"))
    all(isfinite, y) || throw(ArgumentError("QTCF block target values must be finite"))

    (first(z) == last(z) || any(in(nodes(r)), z)) &&
        throw(ArgumentError("interpolation nodes must be distinct"))

    # Reduce both values against the current approximant before fitting the block.
    h = [_new_weight(method, r, z[k], y[k]) for k in eachindex(z)]
    any(isnan, h) && throw(NaNException("Adding QTCF block caused NaN reduced value"))
    all(isfinite, h) || throw(ArgumentError("QTCF block reduced values must be finite"))

    b = _block_line(first(z), first(h), last(z), last(h))
    _qtcf_append_block!(r, z, y, b)
    return r
end

# %% _qtcf_append_block!
# Append accepted data and coefficients, invalidate the cache, and return r.
function _qtcf_append_block!(r::QuadraticThiele, z, y, b)
    S = eltype(r)
    all(isfinite, y) || throw(ArgumentError("QTCF block target values must be finite"))
    any(isnan, b) && throw(NaNException("Adding QTCF block caused NaN weight"))
    all(isfinite, b) || throw(ArgumentError("QTCF block weights must be finite"))

    # Keep scalar storage when possible; otherwise append an explicit block.
    data = getfield(r, :data)
    if data isa Thiele && length(z) == length(b) == 1
        push!(data.nodes, only(z))
        push!(data.values, only(y))
        push!(data.weights, only(b))
    else
        if data isa Thiele
            r.data = data = _QTCFBlocks{S}(r.pairs, r.bpolys, r.apolys, data.nodes, data.values)
        end

        pair = (S(first(z)), S(last(z)))
        push!(data.pairs, pair)
        push!(data.bpolys, S[b...])
        push!(data.apolys, convert(Vector{S}, _factor_poly(pair...)))
        append!(data.nodes, z)
        append!(data.vals, y)
    end

    # Appending a block invalidates the cached numerator and denominator.
    r.polynomials = nothing
    return r
end

# %% _new_weight
# Return the reduced value for the proposed datum using Classic or OneDiv.
function _new_weight(method::Union{Classic,OneDiv}, r::QuadraticThiele, z_new, y_new)
    d = getfield(r, :data)
    d isa Thiele && return _new_weight(method, d.nodes, d.weights, z_new, y_new)
    S = promote_type(eltype(r), typeof(z_new), typeof(y_new))
    if method isa OneDiv
        A, B = one(S), S(y_new)
        for k in eachindex(d.bpolys)
            a = _evaluate_block(d.apolys[k], d.pairs[k], z_new; numerator=true)
            b = _evaluate_block(d.bpolys[k], d.pairs[k], z_new)
            A, B = -b * A + B, a * A
        end
        return B / A
    end

    # Classic reduction applies each stored block directly to the current value.
    u = S(y_new)
    for k in eachindex(d.bpolys)
        pair = d.pairs[k]
        N = _evaluate_block(d.apolys[k], pair, z_new; numerator=true)
        D = u - _evaluate_block(d.bpolys[k], pair, z_new)
        u = _classic_reduced_value(N, D)
    end
    return u
end

# %% _block_line
# Interpolate the two reduced values by a constant or affine polynomial.
function _block_line(z1, h1, z2, h2)
    slope = (h2 - h1) / (z2 - z1)
    midpoint = (z1 + z2) / 2
    midpoint_value = (h1 + h2) / 2
    return [midpoint_value - slope * midpoint, slope]
end

# %% _factor_poly
# Return linear or quadratic numerator coefficients for the supplied node pair.
function _factor_poly(z1, z2)
    S = promote_type(typeof(float(z1)), typeof(float(z2)))
    singleton = z1 == z2
    z1, z2 = S(z1), S(z2)
    return singleton ? S[-z1, 1] : S[z1*z2, -(z1 + z2), 1]
end

# %% [markdown]
# ## 5. Symmetry wrappers and mesh/selection helpers
#
# `Symmetric` and `symmetric` wrap point and value maps for the adapters below.
# `_clean_param`, `_reflection_param`, and `_mesh_params` locate nodes and build
# the continuum mesh.
# `_new_weight` reduces proposed values in both approximation methods.

# %% Symmetric — wrapper type
"""
    Symmetric(object; reflection=nothing)

Wrap a domain or function with a reflection map. Use [`symmetric`](@ref) for
the default segment point reflection or function value reflection.
"""
struct Symmetric{D,R}
    domain::D
    reflection::R
end

# %% Symmetric — constructor
# Wrap an object with its optional reflection map.
Symmetric(d; reflection=nothing) = Symmetric(d, reflection)

# %% symmetric — domain
"""
    symmetric(domain; reflection=nothing)
    symmetric(segment::Segment; reflection=...)
    symmetric(f::Function; reflection=conj)

Wrap a domain's point reflection or a function's value reflection. A segment
defaults to reflection about its midpoint; an explicit map takes precedence.
In `approximate(symmetric(f), symmetric(d))`,
the maps declare `f(point_reflection(z)) == value_reflection(f(z))`.
Both nodes enter the parameter mesh before refinement. Parameter locations use
curve geometry; the interpolation coordinates remain the supplied reflected points.
"""
symmetric(d; kwargs...) = Symmetric(d; kwargs...)

# %% symmetric — function
# Wrap a function with its value-reflection map.
symmetric(f::Function; reflection::Function=conj) = Symmetric(f; reflection)

# %% symmetric — segment
# Wrap a segment with its supplied map or the default midpoint reflection.
function symmetric(s::Segment; kwargs...)
    haskey(kwargs, :reflection) && return invoke(symmetric, Tuple{Any}, s; kwargs...)

    midpoint = point(s, 1//2)
    return invoke(symmetric, Tuple{Any}, s;
        reflection=z -> midpoint - (z - midpoint), kwargs...)
end

# %% _clean_param
# Return a wrapped or clamped curve parameter.
# Specialize on D because curve objects can also be callable.
function _clean_param(d::D, s) where {D}
    iszero(s) && return zero(s)
    isclosed(d) || return min(max(s, zero(s)), oftype(s, length(d)))

    # Closed curves identify both ends of the parameter interval.
    return mod(s, oftype(s, length(d)))
end

# %% _reflection_param
# Locate a paired coordinate without changing the supplied interpolation point.
function _reflection_param(d::D, s, reflection::R; candidate=nothing) where {D,R}
    T, L = typeof(s), oftype(s, length(d))
    target = reflection(point(d, _clean_param(d, s)))
    best, distance = zero(T), Inf
    function consider(t)
        isfinite(t) || return false
        !isclosed(d) && !(zero(T) <= t <= L) && return false
        sr = _clean_param(d, T(t))
        error = abs(point(d, sr)-target)
        if error < distance
            best, distance = sr, error
        end
        return iszero(error)
    end
    consider(zero(T)) && return best
    !isclosed(d) && consider(L) && return best
    candidate !== nothing && consider(candidate) && return best
    inverse = if d isa ComplexRegions.Arc
        mod(ComplexRegions.arg(d.circle, target)-d.start, one(T))/d.delta
    elseif d isa Union{ComplexRegions.Circle,ComplexRegions.Segment,ComplexRegions.Ray}
        ComplexRegions.arg(d, target)
    else
        nothing
    end
    inverse !== nothing && consider(inverse) && return best
    for t in (s, L-s, s+L/T(2), L/T(2)-s)
        consider(t) && return best
    end

    # Improve the parameter using the curve geometry, with no distance cutoff.
    for (k, side) in enumerate(d isa ComplexRegions.AbstractPath ? d : (d,))
        for seed in (zero(T), T(1)/4, T(1)/2, T(3)/4, one(T))
            t = seed
            for _ in 1:4
                isfinite(t) && zero(t) <= t <= one(t) || break
                consider(T(k-1)+t) && return best
                next = t-real((point(side,t)-target)/ComplexRegions.tangent(side,t))
                next == t && break
                t = next
            end
        end
        za, zb = point(side,zero(T))-target, point(side,one(T))-target
        for component in (real,imag)
            fa, fb = component(za), component(zb)
            all(isfinite,(fa,fb)) && !iszero(fa) && !iszero(fb) &&
                signbit(fa) != signbit(fb) || continue
            lo, hi = zero(T), one(T)
            while true
                t = lo+(hi-lo)/2
                (t == lo || t == hi) && break
                consider(T(k-1)+t) && return best
                value = component(point(side,t)-target)
                (!isfinite(value) || iszero(value)) && break
                if signbit(value) == signbit(fa)
                    lo, fa = t, value
                else
                    hi = t
                end
            end
        end
    end
    return best
end

# %% _mesh_params
# Return a sorted mesh retaining endpoints; specialize its arithmetic on T.
function _mesh_params(d::D, params, extra, ::Type{T}) where {D,T}
    L = T(length(d))
    # Retain both endpoints and every distinct cleaned parameter.
    vals = T[zero(T), L]
    for collection in (params, extra), s in collection
        value = _clean_param(d, T(s))
        push!(vals, value)
    end
    return sort!(unique!(vals))
end

# %% _same
# Compare node coordinates exactly.
_same(z, w) = z == w

# Retain the three-argument helper used by existing callers.
_mesh_params(d, extra, ::Type{T}) where {T} = _mesh_params(d, (), extra, T)

# Keep Thiele's interior sample order and skip any sample rounded to the terminal parameter.
function _qtcf_append_interior!(active, path::DiscretizedPath, intervals=eachindex(path.next))
    rows = size(path.params, 1)
    terminal = length(path.path)
    for j in 2:path.nref+1, i in intervals
        index = (j - 1) * rows + i
        path.params[index] < terminal && push!(active, index)
    end
    return active
end

# %% _qtcf_refine!
# Insert a reflected parameter with the same interval split used by Thiele.
function _qtcf_refine!(d::DiscretizedPath, s)
    s = _clean_param(d.path, eltype(d.params)(s))
    for i in eachindex(d.next)
        successor = d.next[i]
        left = d.params[i, 1]
        right = successor == 0 ? length(d.path) : d.params[successor, 1]
        (s == left || s == right) && return nothing
        left < s < right && return _split_interval!(d, i, s)
    end
    throw(KeyError("Indicated new node is not in the discretization"))
end

# %% [markdown]
# ## 6. Continuum and discrete approximation
#
# Separate sampling loops use the same block construction and stopping rules.

# %% _qtcf_float_type
_qtcf_float_type(d::ComplexCurveOrPath) = promote_type(real_type(d), typeof(float(1)))
function _qtcf_float_type(z::AbstractVector)
    F = typeof(float(1))
    for value in z
        F = promote_type(F, typeof(float(real(value))))
    end
    return F
end

# %% _qtcf_allowed
# Resolve the common optional pole restriction on the supplied domain.
function _qtcf_allowed(allowed, geometry, tol)
    allowed === false && return _ -> false
    allowed === :strict || return allowed
    return geometry isa AbstractVector ?
        (z -> minimum(abs(z-t) for t in geometry; init=Inf) > tol) :
        (z -> ComplexRegions.dist(z, geometry) > tol)
end

# %% _qtcf_value
# Evaluate curve targets or retrieve supplied data.
_qtcf_value(s::_QTCFCurveSampler, z) = s.target(z)
_qtcf_value(s::_QTCFDiscreteSampler, i, z) = s.values[i]
_qtcf_value(s::_QTCFCurveSampler, i, z) = _qtcf_value(s, z)

# %% _qtcf_partner
# Locate the supplied partner without changing its coordinate.
function _qtcf_partner(s::_QTCFCurveSampler, id, z)
    candidate = s.parameter_reflection === nothing ? nothing : s.parameter_reflection(id)
    return _reflection_param(s.domain, id, _ -> z; candidate)
end
function _qtcf_partner(s::_QTCFDiscreteSampler, id, z)
    k = findfirst(i -> s.points[i] == z, s.remaining)
    return k === nothing ? nothing : s.remaining[k]
end

# %% _qtcf_block
# Form a chronological block using actual samples and an optional declared value map.
function _qtcf_block(s, id, z, y, reflection, value_reflection; second=nothing)
    isfinite(y) || return nothing
    if reflection === nothing
        second === nothing && return (ids=(id,), nodes=(z,), values=(y,), raw=(y,))
        j, w, v = second
    else
        w = reflection(z)
        w == z && return (ids=(id,), nodes=(z,), values=(y,), raw=(y,))
        j = _qtcf_partner(s, id, w)
        j === nothing && return nothing
        v = _qtcf_value(s, j, w)
    end
    isfinite(v) || return nothing
    target = value_reflection === nothing ? v : value_reflection(z, y, w)
    return (ids=(id,j), nodes=(z,w), values=(y,target), raw=(y,v))
end

# %% _qtcf_start_index
# Infer only an exact midpoint or an unambiguous fixed sample; never move a sample.
function _qtcf_start_index(z, initial_node, reflection)
    isempty(z) && throw(ArgumentError("QTCF requires at least one sample"))
    if initial_node !== nothing
        index = findfirst(==(initial_node), z)
        index === nothing && throw(ArgumentError("initial_node must be present in the grid"))
        return index
    end
    midpoint = if all(t -> t isa Real, z)
        lo, hi = extrema(z)
        lo/2 + hi/2
    else
        realpart = extrema(real, z)
        imagpart = extrema(imag, z)
        complex(realpart[1]/2+realpart[2]/2, imagpart[1]/2+imagpart[2]/2)
    end
    index = findfirst(==(midpoint), z)
    index !== nothing && return index
    if reflection !== nothing
        fixed = unique(t + zero(t) for t in z if t == reflection(t))
        length(fixed) == 1 && return findfirst(==(only(fixed)), z)
    end
    throw(ArgumentError(
        "Supply initial_node when the grid has no exact midpoint or unique fixed point"))
end

# %% _qtcf_fixed_points
# Validate the ordered scalar startup nodes.
function _qtcf_fixed_points(fixed_points, reflection, initial_node)
    fixed_points === nothing && return nothing
    fixed = collect(fixed_points)
    isempty(fixed) && throw(ArgumentError("Supply at least one fixed point"))
    reflection === nothing && throw(ArgumentError("fixed_points requires a reflection map"))
    length(unique(fixed)) == length(fixed) || throw(ArgumentError("Fixed points must be distinct"))
    all(z -> z == reflection(z), fixed) || throw(ArgumentError(
        "Every supplied point must be fixed by reflection"))
    initial_node !== nothing && initial_node != first(fixed) && throw(ArgumentError(
        "initial_node must equal the first supplied fixed point"))
    return fixed
end

# Evaluate a check set into caller-owned error and recurrence buffers.
function _sweep!(errors, rbuf, abuf, bbuf, r::QuadraticThiele, z, y,
    eval_method::ThieleMethod=default_eval_method())
    m = length(z)
    resize!(errors, m)
    resize!(rbuf, m)
    resize!(abuf, m)
    resize!(bbuf, m)
    _evaluate_thiele_samples!(rbuf, r, z, abuf, bbuf, eval_method)
    @inbounds for i in 1:m
        errors[i] = abs(y[i] - rbuf[i])
    end
    return errors
end

# %% approximate — continuum sampling
"""
    approximate(f, curve, QTCF(); fixed_points, initial_node, float_type, tol, max_degree=100, kwargs...)

Use adaptive interior samples, with the terminal parameter excluded from the
greedy search. Startup uses the parameter midpoint unless
initial_node is supplied. With fixed_points, use each supplied fixed point as
one scalar block, in order, before greedy selection. Later symmetric blocks
contain pairs by default; other fixed points are skipped.
Supplied initial and reflected nodes are assumed to lie on the curve.
Refinement decreases
from initial_refinement to refinement, once per completed block. Symmetric blocks
use reflection; value_reflection optionally declares their partner values.
The check set includes selected nodes, and its target scale never decreases.
Every completed block is recorded; return the lowest-error admissible record.
Budgets count blocks.
Input data determine storage; float_type sets the default tolerance.
"""
function approximate(f::Function, d::ComplexCurveOrPath, ::QuadraticThiele;
    float_type::Type{<:AbstractFloat}=_qtcf_float_type(d),
    tol::Real=1000eps(float_type), allowed=true, max_degree::Int=100,
    max_iter::Int=max_degree, stagnation::Int=5, reflection::Union{Nothing,Function}=nothing,
    value_reflection::Union{Nothing,Function}=nothing, initial_node=nothing,
    fixed_points=nothing, pairs_only::Bool=(reflection !== nothing),
    weight_method::ThieleMethod=default_weight_method(),
    eval_method::ThieleMethod=default_eval_method(), refinement::Int=3,
    initial_refinement::Int=15, parameter_reflection::Union{Nothing,Function}=nothing,
    comparison_scale=nothing)
    refinement >= 1 || throw(ArgumentError("refinement must be positive"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be positive"))
    pairs_only && reflection === nothing && throw(ArgumentError("pairs_only requires a reflection map"))
    value_reflection !== nothing && reflection === nothing && throw(ArgumentError(
        "value_reflection requires a reflection map"))
    parameter_reflection !== nothing && reflection === nothing && throw(ArgumentError(
        "parameter_reflection requires a reflection map"))
    fixed = _qtcf_fixed_points(fixed_points, reflection, initial_node)
    fixed !== nothing && (initial_node = first(fixed))
    next_fixed = 2
    T = real_type(d)
    path = DiscretizedPath(d, T[0, length(d)];
        refinement=max(refinement, initial_refinement), maxpoints=max(2max_iter+2, max_iter*refinement))
    midpoint = T(length(d))/2
    fixed_params = fixed === nothing ? T[] :
        [_reflection_param(d, midpoint, _ -> z) for z in fixed]
    fixed === nothing || all(point(d, s) == z for (s, z) in zip(fixed_params, fixed)) ||
        throw(ArgumentError("Every supplied fixed point must lie exactly on the curve"))
    start = initial_node === nothing ? midpoint :
        _reflection_param(d, midpoint, _ -> initial_node)
    z = initial_node === nothing ? point(d, midpoint) : initial_node + zero(point(d, midpoint))
    sampler = _QTCFCurveSampler(d, f, parameter_reflection)
    policy = _qtcf_allowed(allowed, d, tol)
    block = _qtcf_block(sampler, start, z, f(z), reflection, value_reflection)
    block === nothing && error("Initial QTCF block was not available.")
    r = QuadraticThiele([block.nodes...], [block.values...];
        block_sizes=[length(block.nodes)], method=weight_method)
    selected_params = T[block.ids...]
    selected_values = collect(block.raw)
    selected_nodes = Set(w + zero(w) for w in nodes(r))
    reset!(path, _mesh_params(d, (), selected_params, T); refinement=path.nref)

    # Interior samples retain the native matrix-index order.
    τ, num_ref = path.points, path.nref
    rows = size(τ, 1)
    active = Int[]
    sizehint!(active, 2max_iter * refinement)
    fτ = Matrix{eltype([block.raw...])}(undef, size(τ))
    fmax = zero(real_type(eltype(r)))
    S = eltype(r)
    ztest = Vector{typeof(z)}()
    ytest = Vector{eltype(fτ)}()
    rbuf, abuf, bbuf = S[], S[], S[]
    errors = real_type(S)[]
    history = IterationRecord[IterationRecord(r, NaN, missing)]
    reason, best = :max_degree, 1

    terminal = length(d)
    while true
        empty!(active)
        for j in 2:num_ref+1, i in 1:length(path.next)
            index = (j - 1) * rows + i
            0 < path.params[index] < terminal || continue
            sample = τ[index]
            (sample + zero(sample)) in selected_nodes || push!(active, index)
        end
        _sweep!(fτ, f, τ, active)
        empty!(ztest); empty!(ytest)
        append!(ztest, view(τ, active))
        append!(ztest, nodes(r))
        append!(ytest, view(fτ, active))
        append!(ytest, selected_values)
        _sweep!(errors, rbuf, abuf, bbuf, r, ztest, ytest, eval_method)
        history[end].error = maximum(errors)
        fmax = max(fmax, maximum(abs, ytest))
        forced = fixed !== nothing && next_fixed <= length(fixed)
        eligible = pairs_only ?
            [k for k in eachindex(active) if _eligible_thiele_pair(τ[active[k]],
                reflection(τ[active[k]]), selected_nodes)] :
            collect(eachindex(active))
        if forced
            if length(history) >= max_iter
                reason, best = :max_degree, best_acceptable(history, policy)
                break
            end
            s = fixed_params[next_fixed]
            w = fixed[next_fixed] + zero(z)
            block = _qtcf_block(sampler, s, w, f(w), reflection, value_reflection)
        else
            reason, best = quitting_check(history, stagnation, tol, something(comparison_scale, fmax), max_iter, policy)
            reason === :iterating || break
            if isempty(eligible)
                reason, best = :exhausted, best_acceptable(history, policy)
                break
            end
            _, worst = findmax(view(errors, eligible))
            idx = active[eligible[worst]]
            second = nothing
            if reflection === nothing
                next_errors = [k == worst ? -Inf : errors[i]
                    for (k, i) in enumerate(eligible)]
                error, next = findmax(next_errors)
                if error == -Inf
                    reason, best = :node_failure, best_acceptable(history, policy)
                    break
                end
                j = active[eligible[next]]
                second = (path.params[j], τ[j], fτ[j])
            end
            block = _qtcf_block(sampler, path.params[idx], τ[idx], fτ[idx],
                reflection, value_reflection; second)
        end
        if block === nothing
            reason, best = :node_failure, best_acceptable(history, policy)
            break
        end
        try
            add_node!(r, [block.nodes...], [block.values...], weight_method)
        catch error
            error isa Union{ArgumentError,NaNException} || rethrow()
            reason, best = :node_failure, best_acceptable(history, policy)
            break
        end
        forced && (next_fixed += 1)
        push!(history, IterationRecord(r, NaN, missing))
        append!(selected_params, block.ids)
        append!(selected_values, block.raw)
        union!(selected_nodes, (S(w) + zero(S) for w in block.nodes))

        # Rebuild all interval interiors at the next block refinement level.
        try
            num_ref = max(refinement, num_ref - 1)
            reset!(path, _mesh_params(d, (), selected_params, T); refinement=num_ref)
        catch error
            error isa Union{ArgumentError,MaxRefinementException} || rethrow()
            reason, best = :refinement, best_acceptable(history, policy)
            break
        end
    end
    best = best_acceptable(history, policy)
    r = copy(history[best].interpolant)
    stop = ConvergenceStatus(reason, best, history)
    return ContinuumApproximation(f, d, r, policy, path, history, stop)
end

# %% approximate — discrete sampling
"""
    approximate(values, points, QTCF(); fixed_points, initial_node, float_type, tol, max_degree=100, kwargs...)

Fit a fixed grid in its supplied order. Startup uses the exact geometric midpoint,
or the unique fixed sample when the midpoint is absent; otherwise supply
initial_node. With fixed_points, use each supplied fixed point as one scalar
block, in order. Later symmetric blocks contain pairs by default; other fixed
points are skipped. All startup points and partners must be present in the grid.
Input data determine storage; float_type sets the default tolerance.
The check set includes selected nodes. Every completed block is recorded.
Return the current approximant on convergence, otherwise the lowest-error
admissible record. Budgets count blocks.
"""
function approximate(y::AbstractVector{T}, z::AbstractVector{Z}, ::QuadraticThiele;
    float_type::Type{<:AbstractFloat}=_qtcf_float_type(z),
    tol::Real=1000eps(float_type), allowed=true, max_degree::Int=100,
    max_iter::Int=max_degree, stagnation::Int=5, reflection::Union{Nothing,Function}=nothing,
    value_reflection::Union{Nothing,Function}=nothing, initial_node=nothing,
    fixed_points=nothing, pairs_only::Bool=(reflection !== nothing),
    weight_method::ThieleMethod=default_weight_method(),
    eval_method::ThieleMethod=default_eval_method()) where {T<:Number,Z<:Number}
    length(y) == length(z) || throw(DimensionMismatch("node and value counts must agree"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be positive"))
    pairs_only && reflection === nothing && throw(ArgumentError("pairs_only requires a reflection map"))
    value_reflection !== nothing && reflection === nothing && throw(ArgumentError(
        "value_reflection requires a reflection map"))
    fixed = _qtcf_fixed_points(fixed_points, reflection, initial_node)
    fixed !== nothing && (initial_node = first(fixed))
    next_fixed = 2
    fixed === nothing || all(w -> w in z, fixed) ||
        throw(ArgumentError("Every supplied fixed point must be present in the grid"))
    start = _qtcf_start_index(z, initial_node, reflection)
    sampler = _QTCFDiscreteSampler(y, z, collect(eachindex(z)))
    policy = _qtcf_allowed(allowed, z, tol)
    block = _qtcf_block(sampler, start, z[start], y[start], reflection, value_reflection)
    block === nothing && error("Initial QTCF block was not available.")
    r = QuadraticThiele(collect(block.nodes), collect(block.values);
        block_sizes=[length(block.nodes)], method=weight_method)
    filter!(i -> i ∉ block.ids, sampler.remaining)
    history = IterationRecord[IterationRecord(r, NaN, missing)]
    fmax = maximum(abs, y)
    S = promote_type(eltype(r), T, Z)
    rbuf, abuf, bbuf = S[], S[], S[]
    errors = isconcretetype(S) ? real_type(S)[] : Real[]
    reason, best = :max_degree, 1

    while true
        _sweep!(errors, rbuf, abuf, bbuf, r, z, y, eval_method)
        history[end].error = maximum(errors)
        forced = fixed !== nothing && next_fixed <= length(fixed)
        active = pairs_only ? filter(i -> z[i] != reflection(z[i]), sampler.remaining) :
            sampler.remaining
        if forced
            if length(history) >= max_iter
                reason, best = :max_degree, best_acceptable(history, policy)
                break
            end
            id = findfirst(==(fixed[next_fixed]), z)
            block = _qtcf_block(sampler, id, z[id], y[id], reflection, value_reflection)
        else
            reason, best = quitting_check(history, stagnation, tol, fmax, max_iter, policy)
            reason === :iterating || break
            if isempty(active)
                reason, best = :exhausted, best_acceptable(history, policy)
                break
            end
            _, worst = findmax(view(errors, active))
            id = active[worst]
            second = nothing
            if reflection === nothing
                next_errors = [k == worst ? -Inf : errors[i]
                    for (k, i) in enumerate(active)]
                _, next = findmax(next_errors)
                if next_errors[next] == -Inf
                    reason, best = :node_failure, best_acceptable(history, policy)
                    break
                end
                j = active[next]
                second = (j, z[j], y[j])
            end
            block = _qtcf_block(sampler, id, z[id], y[id], reflection, value_reflection; second)
        end
        if block === nothing
            reason, best = :node_failure, best_acceptable(history, policy)
            break
        end
        try
            add_node!(r, collect(block.nodes), collect(block.values), weight_method)
        catch error
            error isa Union{ArgumentError,NaNException} || rethrow()
            reason, best = :node_failure, best_acceptable(history, policy)
            break
        end
        forced && (next_fixed += 1)
        push!(history, IterationRecord(r, NaN, missing))
        filter!(i -> i ∉ block.ids, sampler.remaining)
    end
    r = copy(history[best].interpolant)
    stop = ConvergenceStatus(reason, best, history)
    return DiscreteApproximation(y, z, r, _remaining_mask(r, z), policy, history, stop)
end

# Evaluate supplied functions once before using the fixed-sample adapter.
approximate(f::Function, z::AbstractVector, q::QuadraticThiele; kwargs...) =
    approximate(f.(z), z, q; kwargs...)

# %% approximate — wrapped function default
"""
    approximate(symmetric(f; reflection=conj), symmetric(domain); kwargs...)
    approximate(f, symmetric(domain); kwargs...)

Select QTCF using the wrapped point and value maps. A wrapped function requires an explicit point map or a symmetric domain.
Plain functions evaluate their own partner values. See [`symmetric`](@ref).
"""
function approximate(f::Symmetric{<:Function}, d; kwargs...)
    method, options = haskey(kwargs, :method) ?
        _pop_deprecated_method_kw(kwargs) : (QuadraticThiele(), kwargs)
    return approximate(f, d, method; options...)
end

# %% approximate — wrapped domain default
# Dispatch a wrapped-domain call, using QTCF when no method is supplied.
function approximate(f::Union{Function,AbstractVector}, d::Symmetric; kwargs...)
    method, options = haskey(kwargs, :method) ?
        _pop_deprecated_method_kw(kwargs) : (QuadraticThiele(), kwargs)
    return approximate(f, d, method; options...)
end

# %% approximate — function symmetry adapter
# Validate and unwrap function symmetry before fitting the data.
function approximate(
    f::Symmetric{<:Function}, d, ::QuadraticThiele;
    reflection=nothing, kwargs...
    )
    # Resolve the point map independently of the wrapped value map.
    haskey(kwargs, :value_reflection) && throw(ArgumentError(
        "Supply the value reflection through symmetric(f; reflection=...), not also value_reflection"))
    f.reflection isa Function || throw(ArgumentError(
        "A symmetric function requires a value reflection map"))
    if d isa Symmetric
        reflection === nothing || throw(ArgumentError(
            "Supply the point reflection through symmetric(domain; reflection=...), not also approximate"))
        reflection = d.reflection
        d = d.domain
    end
    reflection isa Function || throw(ArgumentError(
        "A symmetric function requires a point reflection; use symmetric(domain; reflection=...) or the reflection keyword"))

    # Adapt the wrapped value map to the construction callback.
    value_reflection = (z, y, zr) -> f.reflection(y)
    if d isa AbstractVector
        return approximate(f.domain.(d), d, QuadraticThiele(); reflection, value_reflection, kwargs...)
    end
    return approximate(f.domain, d, QuadraticThiele(); reflection, value_reflection, kwargs...)
end
# %% approximate — domain symmetry adapter
# Forward a wrapped domain and its point reflection to QTCF.
function approximate(f::Union{Function,AbstractVector}, d::Symmetric, ::QuadraticThiele; kwargs...)
    f isa AbstractVector && !(d.domain isa AbstractVector) &&
        throw(ArgumentError("sampled values require a vector of domain points"))
    return approximate(f, d.domain, QuadraticThiele();
        reflection=d.reflection, kwargs...)
end

# %% [markdown]
# ## 7. Arithmetic operations
#
# `+`, unary `-`, and `*` match Thiele's scalar operations. Scalar storage delegates
# to Thiele; paired storage changes the first block or scales alternate blocks.

# %% Base.:+
# Return r + s, shifting the first block and stored target values.
function Base.:+(r::QuadraticThiele, s::Number)
    d = getfield(r, :data)
    d isa Thiele && return QuadraticThiele(d + s, nothing)

    S = promote_type(eltype(r), typeof(s))
    b = [S.(p) for p in d.bpolys]
    b[1][1] += s
    return _qtcf_from_blocks(Tuple{S,S}.(d.pairs), b,
        [S.(p) for p in d.apolys], S.(d.nodes), S.(d.vals .+ s))
end

# %% Base.:-
# Return -r through scalar multiplication.
function Base.:-(r::QuadraticThiele)
    return r * -one(eltype(r))
end

# %% Base.:*
# Return r * s, scaling alternate denominator blocks by s and its reciprocal.
function Base.:*(r::QuadraticThiele, s::Number)
    d = getfield(r, :data)
    d isa Thiele && return QuadraticThiele(d * s, nothing)
    S = promote_type(eltype(r), typeof(s))
    if iszero(s)
        z = S(first(d.nodes))
        return QuadraticThiele(S[z], S[0], S[0])
    end

    b = [isodd(j) ? s .* S.(p) : S.(p) ./ s for (j, p) in enumerate(d.bpolys)]
    return _qtcf_from_blocks(Tuple{S,S}.(d.pairs), b,
        [S.(p) for p in d.apolys], S.(d.nodes), S.(s .* d.vals))
end
