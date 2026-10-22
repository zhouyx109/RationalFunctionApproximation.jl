# -*- coding: utf-8 -*-
# %% [markdown]
# # Quadratic Thiele continued fractions
#
# ## 1. Types and accessors
#
# `_QTCFBlocks` stores blocks; `QuadraticThiele` uses Thiele storage for singletons.
# `nodes`, `values`, `degrees`, and `degree` support RFA results and display.
# `convert` changes precision; `copy` saves snapshots; `isreal` checks coefficients.
# `getproperty` builds continuants for degree and reality queries.

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

# %% QuadraticThiele — approximant type
"""
    QuadraticThiele()
    QuadraticThiele(nodes, values, [weights]; block_sizes, method=default_weight_method())

Quadratic Thiele representation, aliased as `QTCF`. Pass the empty selector to
`approximate`, or construct from node-value data. `block_sizes` contains ones
and twos (default: all ones); optional weights are constant or affine coefficient
vectors, one per block. Singleton storage delegates to Thiele. Continuants are lazy.
"""
mutable struct QuadraticThiele{T,S} <: AbstractRationalInterpolant{T,S}
    data::Union{Thiele{T,S},_QTCFBlocks{S}}
    polynomials::Union{Nothing,Tuple{Vector{S},Vector{S}}}
end

# %% Base.convert — numeric precision
"""
    convert(F, r::QuadraticThiele)

Convert the QTCF coefficients and data to floating-point type `F`, preserving
whether the representation is real or complex.
"""
function Base.convert(::Type{F}, r::QuadraticThiele{T,S}) where {F<:AbstractFloat,T,S}
    V = S <: Complex ? Complex{F} : F
    d, p = getfield(r, :data), getfield(r, :polynomials)
    data = d isa Thiele ? convert(F, d) : _QTCFBlocks{V}(
        Tuple{V,V}.(d.pairs), [V.(b) for b in d.bpolys], [V.(a) for a in d.apolys],
        V.(d.nodes), V.(d.vals))
    return QuadraticThiele{F,V}(data, p === nothing ? nothing : (V.(p[1]), V.(p[2])))
end

# %% nodes
"Return the selected QTCF interpolation nodes."
nodes(r::QuadraticThiele) = r.nodes

# %% Base.values
"Return the target values at the QTCF interpolation nodes."
values(r::QuadraticThiele) = r.vals

# %% weights
"Return the constant or affine denominator-block coefficient vectors."
weights(r::QuadraticThiele) = r.bpolys

# %% degrees
"Return the QTCF numerator and denominator degrees before common-factor cancellation."
degrees(r::QuadraticThiele) = map((r.P, r.Q)) do p
    something(findlast(!iszero, p), min(1, length(p))) - 1
end

# %% degree
"Return the QTCF denominator degree before common-factor cancellation."
degree(r::QuadraticThiele) = last(degrees(r))

# %% Base.isreal
# Check whether both continuants have real coefficients.
Base.isreal(r::QuadraticThiele) = eltype(r) <: Real || (all(isreal, r.P) && all(isreal, r.Q))

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
# Return stored nodes/values or build the continuants used by degrees.
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
# `evaluate` and `evaluate!` use the shared OneDiv/Classic strategies. `_evaluate_numden` and
# `_evaluate_classic` supply their block recurrences; `_evaluate_block` evaluates
# one block, and `_qtcf_ratio` handles exceptional divisions.
# `_padd` and `_pmul` build continuants in `getproperty`.
# `_evaluate_block_derivs` and `_evaluate_numden_derivs` supply polynomial
# derivatives for `residues`; `_derivative!` supplies general quotient derivatives.

# %% evaluate — scalar
"""
    evaluate(r::QuadraticThiele, z, method=default_eval_method())

Evaluate one point with `OneDiv()` or backward `Classic()` evaluation.
OneDiv retries Classic at a zero denominator; scalar evaluation guards underflow.
At a selected common zero, evaluate after canceling its simple common factor.
"""
function evaluate(r::QuadraticThiele, z::Number, method::ThieleMethod=default_eval_method())
    d = getfield(r, :data)
    d isa Thiele && !isempty(d.weights) && !isfinite(z) &&
        return convert(promote_type(eltype(r), typeof(z)), evaluate(d, z, method))

    # Resolve nonfinite inputs before entering the block recurrence.
    if isinf(z)
        P, Q = r.P, r.Q
        m, n = degrees(r)
        return m > n ? oftype(last(P), Inf) :
            m < n ? zero(last(P)) : _qtcf_ratio(last(P), last(Q))
    elseif isnan(z)
        return oftype(z, NaN)
    end

    S = promote_type(eltype(r), typeof(z))
    value = d isa Thiele && !isempty(d.weights) ?
        convert(S, evaluate(d, z, method)) : _evaluate(r, z, method)

    # Cancel a selected-node 0/0 using the first derivatives of the continuants.
    if !isfinite(value) && any(==(z), nodes(r))
        A, B = ntuple(_ -> Vector{S}(undef, 2), 2)
        _derivative!(z, r, 1, A, B, A; quotient=false)
        if iszero(A[1]) && iszero(B[1]) && all(isfinite, (A[2], B[2])) &&
                !(iszero(A[2]) && iszero(B[2]))
            return _qtcf_ratio(A[2], B[2])
        end
    end
    return value
end

# %% evaluate — array
# Evaluate an array of points into a promoted output array.
function evaluate(r::QuadraticThiele, z::AbstractArray{<:Number}, method::ThieleMethod=default_eval_method())
    t = similar(z, promote_type(eltype(r), eltype(z)))
    return evaluate!(t, r, z, method)
end

# %% evaluate! — allocate scratch
"""
    evaluate!(t, r::QuadraticThiele, z, [a, b], method=default_eval_method())

Write evaluations into `t`. Scratch arrays must match axes and must not overlap
one another, input, or output. Input/output aliasing is supported; arrays do not rescale.
"""
function evaluate!(
    t::AbstractArray, r::QuadraticThiele, z::AbstractArray{<:Number},
    method::ThieleMethod=default_eval_method()
    )
    return evaluate!(t, r, z, similar(t), similar(t), method)
end

# %% evaluate! — reuse scratch
# Fill t using caller-supplied scratch arrays a and b.
function evaluate!(
    t::AbstractArray, r::QuadraticThiele, z::AbstractArray{<:Number},
    a::AbstractArray, b::AbstractArray, method::ThieleMethod=default_eval_method()
    )
    # Validate scratch storage and preserve an aliased input.
    axes(t) == axes(z) == axes(a) == axes(b) ||
        throw(DimensionMismatch("evaluation input, output, and scratch arrays must have matching axes"))
    (Base.mightalias(a, b) || Base.mightalias(a, t) || Base.mightalias(b, t) ||
        Base.mightalias(a, z) || Base.mightalias(b, z)) &&
        throw(ArgumentError("evaluation scratch arrays must not overlap each other, the input, or the output"))
    source = Base.mightalias(t, z) ? copy(z) : z

    # Singleton storage uses the scalar TCF kernels.
    d, S = getfield(r, :data), promote_type(eltype(r), eltype(source))
    if d isa Thiele
        isempty(d.weights) && return map!(x -> evaluate(r, x, method), t, source)
        method isa OneDiv || return map!(x -> evaluate(r, x, method), t, source)
        evaluate!(t, d, source, a, b)
        for i in eachindex(t, source)
            (isnan(source[i]) || !isfinite(t[i])) && (t[i] = evaluate(r, source[i], method))
        end
        return t
    end

    method isa OneDiv || return map!(x -> evaluate(r, x, method), t, source)
    m = length(d.bpolys)

    # Seed the pair recurrence from the final block.
    @inbounds for i in eachindex(a, b, source)
        a[i] = S(_evaluate_block(d.bpolys[m], d.pairs[m], source[i]))
        b[i] = m == 1 ? one(S) :
            S(_evaluate_block(d.apolys[m-1], d.pairs[m-1], source[i]; numerator=true))
    end

    # Propagate through the interior blocks without division.
    @inbounds for k in m-1:-1:2
        block, pair = d.bpolys[k], d.pairs[k]
        affine = length(block) == 2
        half = oftype(float(real(pair[1])), 0.5)
        midpoint = half * pair[1] + half * pair[2]
        slope = affine ? block[2] : zero(eltype(block))
        value = affine ? block[1] + slope * midpoint : block[1]
        numerator, (z1, z2) = d.apolys[k-1], d.pairs[k-1]
        quadratic, leading = length(numerator) == 3, last(numerator)
        for i in eachindex(a, b, source)
            A = a[i]
            bj = affine ? value + slope * (source[i] - midpoint) : value
            aj = leading * (source[i] - z1)
            quadratic && (aj *= source[i] - z2)
            a[i], b[i] = bj * A + b[i], A * aj
        end
    end

    # Apply the first block and retry Classic at a zero denominator.
    @inbounds for i in eachindex(t, a, b, source)
        A = a[i]
        if m > 1
            bj = _evaluate_block(d.bpolys[1], d.pairs[1], source[i])
            a[i], b[i] = bj * A + b[i], A
        end
        t[i] = !isfinite(source[i]) ? evaluate(r, source[i], method) :
            iszero(b[i]) ? evaluate(r, source[i], method) : a[i] / b[i]
    end
    return t
end

# %% _evaluate_classic
# Evaluate by the backward continued-fraction recurrence.
function _evaluate_classic(r::QuadraticThiele, z::Number)
    d = getfield(r, :data)
    d isa Thiele && return isempty(d.weights) ?
        zero(promote_type(eltype(r), typeof(z))) : _evaluate_classic(d, z)
    m = length(d.bpolys)
    u = _evaluate_block(d.bpolys[m], d.pairs[m], z)
    for k in m-1:-1:1
        u = _evaluate_block(d.bpolys[k], d.pairs[k], z) +
            _qtcf_ratio(_evaluate_block(d.apolys[k], d.pairs[k], z; numerator=true), u)
    end
    return u
end

# %% _evaluate_numden
# Return the backward numerator/denominator pair with underflow scaling.
function _evaluate_numden(r::QuadraticThiele, z::Number; underflow_guard=true)
    @assert isfinite(z)
    d, S = getfield(r, :data), promote_type(eltype(r), typeof(z))
    if d isa Thiele
        isempty(d.weights) && return zero(S), one(S)
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

# %% _evaluate_block
# Evaluate a constant/affine denominator or a factored numerator block at z.
@inline function _evaluate_block(p, pair, z; numerator=false)
    if numerator
        return length(p) == 2 ? p[2] * (z - pair[1]) :
            p[3] * (z - pair[1]) * (z - pair[2])
    end
    length(p) == 1 && return only(p)

    # Center affine evaluation on the interpolation pair.
    half = oftype(float(real(pair[1])), 0.5)
    midpoint = half * pair[1] + half * pair[2]
    midpoint_value = p[1] + p[2] * midpoint
    return midpoint_value + p[2] * (z - midpoint)
end

# %% _qtcf_ratio
# Return numer / denom with explicit NaN, Inf, and zero-denominator handling.
function _qtcf_ratio(numer, denom)
    n, d = promote(numer, denom)
    (isnan(n) || isnan(d)) && return oftype(d, NaN)
    iszero(d) && return oftype(d, iszero(n) ? NaN : Inf)
    isinf(d) && return oftype(d, isinf(n) ? NaN : 0)
    return n / d
end

# %% _padd
# Return the trimmed sum of two polynomial coefficient vectors.
function _padd(p, q)
    n = max(length(p), length(q))
    r = zeros(promote_type(eltype(p), eltype(q)), n)
    @views r[1:length(p)] .+= p
    @views r[1:length(q)] .+= q
    return resize!(r, something(findlast(x -> !(abs(x) <= 0), r), min(1, length(r))))
end

# %% _pmul
# Return the trimmed product of two polynomial coefficient vectors.
function _pmul(p, q)
    r = zeros(promote_type(eltype(p), eltype(q)), length(p) + length(q) - 1)
    @inbounds for i in eachindex(p), j in eachindex(q)
        r[i+j-1] += p[i] * q[j]
    end
    return resize!(r, something(findlast(x -> !(abs(x) <= 0), r), min(1, length(r))))
end

# %% _evaluate_block_derivs
# Return a block value and its first two derivatives.
@inline function _evaluate_block_derivs(p, pair, z; numerator=false)
    S = promote_type(ComplexF64, eltype(p), typeof(z))
    if !numerator
        return S(_evaluate_block(p, pair, z)), length(p) > 1 ? S(p[2]) : zero(S), zero(S)
    elseif length(p) == 2
        c, z1 = S(p[2]), S(pair[1])
        return c * (z - z1), c, zero(S)
    end
    c, z1, z2 = S(p[3]), S(pair[1]), S(pair[2])
    return c * (z - z1) * (z - z2), c * ((z - z1) + (z - z2)), 2c
end

# %% _derivative! — workspace and quotient recovery
# Fill numerator/denominator derivatives and recover vals with one reciprocal.
function _derivative!(
    z, r::QuadraticThiele, order::Integer,
    A::AbstractVector, B::AbstractVector, vals::AbstractVector;
    underflow_guard=true, quotient=true
    )
    order >= 0 || throw(ArgumentError("derivative order must be nonnegative"))
    n = order + 1
    min(length(A), length(B), length(vals)) >= n ||
        throw(DimensionMismatch("derivative workspaces must contain at least order + 1 entries"))
    Base.require_one_based_indexing(A, B, vals)
    (Base.mightalias(A, B) || quotient &&
        (Base.mightalias(A, vals) || Base.mightalias(B, vals))) &&
        throw(ArgumentError("derivative workspaces must not overlap"))
    d, S = getfield(r, :data), promote_type(eltype(r), typeof(z))
    if d isa Thiele
        _derivative!(z, d.weights, d.nodes, order, A, B, vals;
            underflow_guard, quotient=false)
    else
        # Seed the final block, then differentiate backward in descending order.
        m = length(d.bpolys)
        fill!(view(A, 1:n), 0); fill!(view(B, 1:n), 0)
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
            if k > 1 && underflow_guard && any(guard, view(A, 1:n))
                view(A, 1:n) .*= 1e50
                view(B, 1:n) .*= 1e50
            end
        end
    end
    quotient || return A, B

    # Recover all quotient derivatives with one reciprocal.
    inv_denom = inv(B[1])
    vals[1] = A[1] * inv_denom
    for q in 1:order
        correction = zero(eltype(vals))
        for j in 1:q
            correction += binomial(q, j) * B[j+1] * vals[q-j+1]
        end
        vals[q+1] = (A[q+1] - correction) * inv_denom
    end
    return vals
end

# %% _evaluate_numden_derivs
# Return unscaled P, Q, P', Q', reusing scratch when supplied.
function _evaluate_numden_derivs(r::QuadraticThiele, z::Number, space=nothing)
    S = promote_type(eltype(r), typeof(z))
    A, B = space === nothing ?
        ntuple(_ -> Vector{promote_type(ComplexF64, S)}(undef, 2), 2) : space
    _derivative!(z, r, 1, A, B, A; underflow_guard=false, quotient=false)
    return S(A[1]), S(B[1]), S(A[2]), S(B[2])
end

# %% derivative — single order
"""
    derivative(r::QuadraticThiele, order=1)
    derivative(r::QuadraticThiele, orders::AbstractVector)

Return an evaluator for the requested nonnegative derivative order(s).
Scalar storage uses Thiele's implementation; block derivatives use one reciprocal.
"""
function derivative(r::QuadraticThiele, order::Integer=1)
    f = derivative(r, [order])
    return z -> only(f(z))
end

# %% derivative — requested orders
# Return requested derivatives, using separate workspaces for each evaluation.
function derivative(r::QuadraticThiele, orders::AbstractVector{<:Integer})
    isempty(orders) && throw(ArgumentError("at least one derivative order is required"))
    minimum(orders) >= 0 || throw(ArgumentError("derivative orders must be nonnegative"))
    d = getfield(r, :data)
    d isa Thiele && return derivative(d, orders)
    order, index = maximum(orders), 1 .+ orders
    return function(z)
        S = promote_type(eltype(r), typeof(z))
        A, B, vals = ntuple(_ -> Vector{S}(undef, order + 1), 3)
        _derivative!(z, r, order, A, B, vals)
        return vals[index]
    end
end

# %% [markdown]
# ## 3. Poles, roots, residues, and partial fractions
#
# `roots` assembles an unreduced numerator or denominator pencil; `poles` cancels
# exact selected-node common factors from the denominator eigenvalues.
# `residues` reuses RFA's residue loop and contour fallback with `_evaluate_numden_derivs`.
# `convert(PartialFractions, r)` combines these with the polynomial part.

# %% poles
"""
    poles(r::QuadraticThiele)

Return finite poles after canceling exact common factors at earlier block nodes.
The unreduced eigensolve delegates to Thiele for scalar storage. Use
`roots(r, 2)` for the finite roots of the unreduced denominator pencil.
"""
poles(r::QuadraticThiele) = first(_qtcf_pole_data(r))

# %% roots
"""
    roots(r::QuadraticThiele, first_block=1)

Return finite roots of the unreduced numerator (`first_block=1`) or denominator
(`first_block=2`) from its block matrix pencil, before common-factor cancellation.
"""
function roots(r::QuadraticThiele, first_block::Integer=1)
    data = getfield(r, :data)
    data isa Thiele && return first_block == 2 ? poles(data) : roots(data)
    bpolys, apolys, pairs = r.bpolys, r.apolys, r.pairs

    # Assemble the constant and linear parts of the block matrix pencil.
    S = promote_type(eltype(r), eltype(r.nodes))
    n = max(length(bpolys) - first_block + 1, 0)
    C, D = zeros(S, n, n), zeros(S, n, n)
    for j in 1:n
        k = j + first_block - 1
        b = bpolys[k]
        C[j, j] += S(b[1])
        length(b) >= 2 && (D[j, j] += S(b[2]))

        # Split each numerator block across the two adjacent off-diagonals.
        if j >= 2
            a = apolys[k-1]
            if length(a) <= 2
                C[j-1, j] += S(a[1])
                length(a) == 2 && (D[j-1, j] += S(a[2]))
                C[j, j-1] += -one(S)
            else
                z1, z2 = S.(pairs[k-1])
                C[j-1, j] += -z1
                D[j-1, j] += one(S)
                C[j, j-1] += z2
                D[j, j-1] += -one(S)
            end
        end
    end

    # Extract finite eigenvalues, with a complex Schur fallback.
    isempty(C) && return eltype(r.nodes)[]

    # Equilibrate both matrices together without changing the pencil roots.
    for _ in 1:40
        changed = false
        for dim in 1:2, j in 1:n
            c = dim == 1 ? view(C, j, :) : view(C, :, j)
            d = dim == 1 ? view(D, j, :) : view(D, :, j)
            scale = max(maximum(abs, c), maximum(abs, d))
            (iszero(scale) || !isfinite(scale)) && continue
            exponent = round(Int, -log2(scale) / 2)
            iszero(exponent) && continue
            factor = ldexp(one(scale), exponent)
            c .*= factor
            d .*= factor
            changed = true
        end
        changed || break
    end
    return try
        filter(isfinite, eigvals(-C, D))
    catch
        _, _, _, _, α, β = schur(complex(-C), complex(D))
        filter(isfinite, α ./ β)
    end
end

# %% _qtcf_pole_data
# Confirm selected-node common factors exactly, then remove one copy of each.
function _qtcf_pole_data(r::QuadraticThiele)
    ζ, common = roots(r, 2), eltype(r)[]
    orders = Int[]
    bpolys, apolys, points = r.bpolys, r.apolys, r.pairs
    m = length(bpolys)
    (m < 2 || isempty(ζ)) && return ζ, common

    # Stored floating-point numbers are exact rationals. A finite-field screen
    # rejects ordinary tails cheaply; only its zeros require a rational recurrence.
    # 256^2 = -1 mod 65537 also permits complex coefficients in the same screen.
    E = eltype(r) <: Real ? Rational{BigInt} : Complex{Rational{BigInt}}
    b, pairs = [E.(p) for p in bpolys], [E.(p) for p in points]
    leading = E[last(a) for a in apolys]
    prime = 65537
    image_real(x) = mod(Int(mod(numerator(x), prime)) *
        invmod(Int(mod(denominator(x), prime)), prime), prime)
    image(x) = mod(image_real(real(x)) + 256image_real(imag(x)), prime)
    mb, mpairs, ml = [image.(p) for p in b], [image.(p) for p in pairs], image.(leading)

    # At a node of block j, a_j vanishes and the lower-tail numerator is zero
    # exactly when P_m and Q_m share its simple factor (the determinant identity).
    function tail(z, j, blocks, points, coefficients, reduce)
        A, B = reduce(evalpoly(z, last(blocks))), one(z)
        for k in m-1:-1:j+1
            a = coefficients[k] * (z - points[k][1])
            length(apolys[k]) == 3 && (a *= z - points[k][2])
            A, B = reduce(reduce(evalpoly(z, blocks[k])) * A + reduce(a) * B), A
        end
        return A
    end
    for j in 1:m-1, z in (length(apolys[j]) == 2 ?
            (first(points[j]),) : unique(points[j]))
        iszero(tail(image(E(z)), j, mb, mpairs, ml, x -> mod(x, prime))) || continue
        iszero(tail(E(z), j, b, pairs, leading, identity)) || continue
        push!(common, z)

        # Exact Taylor coefficients identify any true pole left at this node.
        # This also avoids returning a split approximation to a known multiple root.
        zz, Q, Qold = E(z), E[1], E[0]
        for k in 2:m
            bj = length(b[k]) == 1 ? b[k] : E[evalpoly(zz, b[k]), b[k][2]]
            z1, z2 = pairs[k-1]
            aj = length(apolys[k-1]) == 2 ? E[zz-z1, 1] :
                E[(zz-z1)*(zz-z2), (zz-z1)+(zz-z2), 1]
            Q, Qold = _padd(_pmul(bj, Q), _pmul(leading[k-1] .* aj, Qold)), Q
        end
        push!(orders, something(findfirst(!iszero, Q), length(Q) + 1) - 1)
    end
    isempty(common) && return ζ, common

    # The exact factor test determines the multiplicity removed. In particular,
    # a close genuine pole is retained; proximity alone never cancels a pole.
    for z in common
        isempty(ζ) && break
        deleteat!(ζ, argmin(abs.(ζ .- z)))
    end
    available = collect(eachindex(ζ))
    for (z, order) in zip(common, orders), _ in 1:order-1
        isempty(available) && break
        j = argmin(abs.(ζ[available] .- z))
        ζ[available[j]] = z
        deleteat!(available, j)
    end
    return ζ, common
end

# %% residues
"""
    residues(r::QuadraticThiele)

Return `(poles, residues)`, using the numerator/denominator derivative ratio
with contour integration after cancellation and as a fallback for nonfinite results.
"""
function residues(r::QuadraticThiele)
    ζ, common = _qtcf_pole_data(r)
    S = promote_type(ComplexF64, eltype(r), eltype(ζ))
    space = isempty(ζ) ? nothing : ntuple(_ -> Vector{S}(undef, 2), 2)
    # Splitting a multiple unreduced denominator root can make P/Q' unreliable
    # even when finite. Contours avoid this after any confirmed cancellation.
    numden(r, z) = !isempty(common) ?
        (zero(S), zero(S), zero(S), zero(S)) : _evaluate_numden_derivs(r, z, space)
    radius(z) = minimum((abs(z - t) / 4 for t in ζ if t != z);
        init=max(one(abs(z)), abs(z)) * 0.25)
    isempty(common) && return residues(r, ζ; check_finite=true, radius, numden)
    return setprecision(BigFloat, max(256, precision(real_type(eltype(r))))) do
        # Wider arithmetic prevents cancellation in evaluations on small contours.
        residues(convert(BigFloat, r), ζ; check_finite=true, radius, numden)
    end
end

# %% Base.convert — partial fractions
"""
    convert(PartialFractions, r::QuadraticThiele)

Convert a QTCF with simple computed poles to polynomial-plus-partial-fraction form.
"""
function Base.convert(::Type{PartialFractions}, r::QuadraticThiele)
    # Divide the numerator by the denominator to obtain the polynomial part.
    remainder, q = copy(r.P), copy(r.Q)
    resize!(remainder, something(findlast(!iszero, remainder), 1))
    resize!(q, something(findlast(!iszero, q), 1))
    all(iszero, q) && throw(DivideError())
    S = promote_type(eltype(remainder), eltype(q))
    dp, dq = length(remainder) - 1, length(q) - 1
    quotient = dp < dq ? S[0] : zeros(S, dp - dq + 1)
    @inbounds for shift in (dp-dq):-1:0
        coefficient = remainder[dq + shift + 1] / q[dq + 1]
        quotient[shift + 1] = coefficient
        for k in 0:dq
            remainder[k + shift + 1] -= coefficient * q[k + 1]
        end
    end
    resize!(quotient, something(findlast(!iszero, quotient), 1))

    # Represent the polynomial part in the shared Arnoldi basis.
    S = promote_type(ComplexF64, eltype(quotient), eltype(r.nodes))
    polynomial = if all(iszero, quotient)
        ArnoldiPolynomial(S[0], ArnoldiBasis(S[0], 0; max_degree=0))
    else
        d = length(quotient) - 1
        center = isempty(r.nodes) ? zero(S) : sum(S, r.nodes) / length(r.nodes)
        radius = max(1.0, maximum(z -> Float64(abs(S(z) - center)), r.nodes; init=0.0))
        sample = d == 0 ? S[center] :
            S[center + radius * cis(2pi * k / (d + 1)) for k in 0:d]
        ArnoldiBasis(sample, d; max_degree=d) \ [evalpoly(z, quotient) for z in sample]
    end

    ζ, residue = residues(r)
    any(ζ[i] == ζ[j] for i in eachindex(ζ) for j in firstindex(ζ):i-1) &&
        throw(ArgumentError("PartialFractions conversion requires simple poles"))
    return PartialFractions(polynomial, ζ, residue)
end

# %% [markdown]
# ## 4. Data construction, node updates, and blocks
#
# Constructors and `add_node!` reuse `_new_weight` for reduced values and
# `_qtcf_append_block!` for storage; `_is_selected_node` also supports RFA's `rewind`.
# Both `approximate` methods use `_qtcf_from_blocks` to select storage,
# `_block_line` to fit denominator coefficients, and `_factor_poly` for numerators.
# The same block helpers serve paired `add_node!` updates.

# %% QuadraticThiele — requested precision
"""
    QuadraticThiele{T}(nodes, values, [weights]; kwargs...)

Construct a QTCF at floating-point precision `T`, checking converted data.
"""
function QuadraticThiele{T}(args...; kwargs...) where {T<:AbstractFloat}
    original = QuadraticThiele(args...; kwargs...)
    r = real_type(eltype(original)) === T ? original : convert(T, original)

    # Validate the converted data, since conversion can merge distinct nodes.
    all(isfinite, nodes(r)) || throw(ArgumentError("interpolation nodes must be finite"))
    all(isfinite, values(r)) || throw(ArgumentError("QTCF block target values must be finite"))
    d = getfield(r, :data)
    (d isa Thiele ? all(isfinite, d.weights) : all(b -> all(isfinite, b), d.bpolys)) ||
        throw(ArgumentError("QTCF block weights must be finite"))

    any(_same(r.nodes[i], r.nodes[j]) for i in eachindex(r.nodes)
        for j in firstindex(r.nodes):i-1) &&
        throw(ArgumentError("interpolation nodes must remain distinct after conversion"))
    return r
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
        S = promote_type(map(T -> T <: Number ? float(T) : Float64, types)...)
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

    # Reject nonfinite or repeated nodes before constructing any blocks.
    all(isfinite, z) || throw(ArgumentError("interpolation nodes must be finite"))
    bound = 100eps(Float64) * maximum(abs, z; init=1)
    any(abs(real(z[i]) - real(z[j])) <= bound && abs(imag(z[i]) - imag(z[j])) <= bound &&
        _same(z[i], z[j]) for i in eachindex(z) for j in firstindex(z):i-1) &&
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
# Wrap owned arrays in scalar storage when possible, otherwise in block storage.
function _qtcf_from_blocks(pairs, bpolys, apolys, z::Vector{S}, y) where {S}
    scalar = !isempty(z) && length(z) == length(pairs) == length(bpolys) == length(apolys) &&
        all(eachindex(z)) do j
            pairs[j][1] == pairs[j][2] == z[j] && length(bpolys[j]) == 1 &&
                isfinite(only(bpolys[j])) && length(apolys[j]) == 2 &&
                apolys[j][1] == -z[j] && apolys[j][2] == 1
        end

    data = scalar ? Thiele(z, y, S[only(b) for b in bpolys]) :
        _QTCFBlocks{S}(pairs, bpolys, apolys, z, y)
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
    _is_selected_node(r, z) && throw(ArgumentError("interpolation nodes must be distinct"))

    h = _new_weight(method, r, z, y)
    isnan(h) && throw(NaNException("Adding QTCF block caused NaN reduced value"))
    isfinite(h) || throw(ArgumentError("QTCF block reduced values must be finite"))
    return _qtcf_append_block!(r, (z,), (y,), (h,))
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

    # Reject nodes already selected or coinciding within the pair.
    any(t -> any(w -> _same(t, w), r.nodes), z) &&
        throw(ArgumentError("interpolation nodes must be distinct"))
    _same(first(z), last(z)) &&
        throw(ArgumentError("a two-point block requires distinct nodes"))

    # Reduce both values against the current approximant before fitting the block.
    h = method isa OneDiv && getfield(r, :data) isa _QTCFBlocks ?
        convert(Vector{eltype(r)}, _reduced_values(z, y, r.bpolys, r.apolys, r.pairs)) :
        [_new_weight(method, r, z[k], y[k]) for k in eachindex(z)]
    any(isnan, h) && throw(NaNException("Adding QTCF block caused NaN reduced value"))
    all(isfinite, h) || throw(ArgumentError("QTCF block reduced values must be finite"))

    b = _block_line(first(z), first(h), last(z), last(h); symmetric=false)
    return _qtcf_append_block!(r, z, y, b)
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

# %% _is_selected_node
# Test whether the proposed node matches a stored node.
function _is_selected_node(r::QuadraticThiele, z)
    z = convert(eltype(r), z)
    scale = max(1, abs(real(z)) + abs(imag(z)))
    return any(nodes(r)) do w
        bound = 200eps(Float64) * max(scale, abs(real(w)) + abs(imag(w)))
        abs(real(z) - real(w)) <= bound && abs(imag(z) - imag(w)) <= bound && _same(z, w)
    end
end

# %% _new_weight
# Return the reduced value for the proposed datum using Classic or OneDiv.
function _new_weight(method::Union{Classic,OneDiv}, r::QuadraticThiele, z_new, y_new)
    d = getfield(r, :data)
    d isa Thiele && return _new_weight(method, d.nodes, d.weights, z_new, y_new)
    S = promote_type(eltype(r), typeof(z_new), typeof(y_new))
    method isa OneDiv && return S(only(_reduced_values((z_new,), (y_new,), r.bpolys, r.apolys, r.pairs)))

    # Classic reduction applies each stored block directly to the current value.
    u = S(y_new)
    for k in eachindex(r.bpolys)
        pair = r.pairs[k]
        a = _evaluate_block(r.apolys[k], pair, z_new; numerator=true)
        b = _evaluate_block(r.bpolys[k], pair, z_new)
        u = _qtcf_ratio(a, u - b)
    end
    return u
end

# %% _block_line
# Fit constant/affine weights; optionally preserve symmetry and retry overflow.
function _block_line(z1, h1, z2, h2; symmetric=true)
    (!isfinite(h1) || !isfinite(h2)) && return [oftype(h1, NaN)]
    _same(z1, z2) && return [h1]

    # Fit the pair while retaining equal or conjugate reduced values.
    b = if symmetric && abs(h2 - h1) <=
            100eps(Float64) * max(abs(h1), abs(h2))
        S = promote_type(typeof(h1), typeof(h2))
        S[(S(h1) + S(h2)) / S(2)]
    elseif symmetric && !iszero(imag(z1)) && _same(z2, conj(z1)) &&
            abs(h2 - conj(h1)) <= 100eps(Float64) * max(abs(h1), abs(h2))
        # Conjugate data determine a line with real coefficients.
        R = promote_type(typeof(real(z1)), typeof(real(h1)))
        γ = R(imag(h1) / imag(z1))
        β = R(real(h1) - γ * real(z1))
        R[β, γ]
    else
        slope = (h2 - h1) / (z2 - z1)
        midpoint = (z1 + z2) / 2
        midpoint_value = (h1 + h2) / 2
        [midpoint_value - slope * midpoint, slope]
    end

    all(isfinite, b) && return b
    S = eltype(b)
    all(isfinite, (z1, h1, z2, h2)) || return S[NaN]

    # Retry finite data in wider arithmetic when the block calculation overflows.
    wide = setprecision(BigFloat, max(256, precision(BigFloat))) do
        if length(b) == 1
            [(big(h1) + big(h2)) / 2]
        else
            wide_slope = (big(h2) - big(h1)) / (big(z2) - big(z1))
            [(big(h1) + big(h2) - wide_slope * (big(z1) + big(z2))) / 2, wide_slope]
        end
    end

    # Wider arithmetic helps only if the coefficients fit the original type.
    result = S.(wide)
    return all(isfinite, result) ? result : S[NaN]
end

# %% _factor_poly
# Return linear or quadratic numerator coefficients for the supplied node pair.
function _factor_poly(z1, z2)
    S = promote_type(ComplexF64, typeof(z1), typeof(z2))
    singleton = _same(z1, z2)
    z1, z2 = S(z1), S(z2)
    return singleton ? S[-z1, 1] : S[z1*z2, -(z1 + z2), 1]
end

# %% [markdown]
# ## 5. Symmetry wrappers and mesh/selection helpers
#
# `Symmetric` and `symmetric` wrap point and value maps for the adapters below.
# `_clean_param`, `_reflection_param`, and `_mesh_params` locate and refine
# continuum nodes; `_reflected_continuum_samples` completes its test samples.
# `_same`, `_node_hash_key`, and `_lookup_node_index` match nodes in both fits.
# `_reduced_values` uses `_projective_reduce!` to reduce proposed continuum
# values; discrete fitting also updates its full reduced grid with that helper.
# `_qtcf_converged` applies the error and pole tests in both approximation loops.

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
Both nodes enter the parameter mesh before refinement. Reflected parameters use
analytic inverses, verified parameter correspondences, or the curve's locator.
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
    isclosed(d) || return min(max(s, zero(s)), oftype(s, length(d)))

    # Closed curves identify both ends of the parameter interval.
    L = convert(typeof(s), length(d))
    t = mod(s, L)
    tol = 100eps(typeof(t)) * max(one(t), abs(L))
    return (t <= tol || abs(t - L) <= tol) ? zero(t) : t
end

# %% _reflection_param
# Return the curve parameter of the supplied reflected point.
function _reflection_param(d::D, s, reflection::R) where {D,R}
    T, L = typeof(s), oftype(s, length(d))
    target = reflection(point(d, _clean_param(d, s)))

    # Check a located parameter against the reflected point and rounding tolerance.
    function matches(t)
        z = point(d, t)
        _same(z, target) && return true
        h = T(4) * eps(T) * max(one(T), abs(t))
        bound = max(abs(point(d, _clean_param(d, t-h)) - z),
            abs(point(d, _clean_param(d, t+h)) - z))
        return isfinite(bound) && abs(z - target) <= bound
    end

    # Try exact endpoints and a curve-provided inverse first.
    _same(point(d, zero(T)), target) && return zero(T)
    !isclosed(d) && _same(point(d, L), target) && return L
    if applicable(ComplexRegions.arg, d, target)
        sr = _clean_param(d, T(ComplexRegions.arg(d, target)))
        matches(sr) && return sr
    end

    # These common parameter correspondences are accepted only after checking the point map.
    for t in (s, L - s, s + L / T(2), L / T(2) - s)
        !isclosed(d) && !(zero(T) <= t <= L) && continue
        sr = _clean_param(d, t)
        matches(sr) && return sr
    end

    # Fall back to locating the target on each side of the path.
    for (k, side) in enumerate(d isa ComplexRegions.AbstractPath ? d : (d,))
        t = applicable(ComplexRegions.arg, side, target) ? ComplexRegions.arg(side, target) :
            last(ComplexRegions.finddist(target, side))

        # Polish the locator in parameter space, retaining the pointwise acceptance check.
        for _ in 1:4
            isfinite(t) && zero(t) <= t <= one(t) || break
            sr = _clean_param(d, T(k - 1) + T(t))
            matches(sr) && return sr
            t -= real((point(side, t) - target) / ComplexRegions.tangent(side, t))
        end

        # A bracketed coordinate root supplies another candidate to verify geometrically.
        za, zb = point(side, zero(T)) - target, point(side, one(T)) - target
        for component in (real, imag)
            fa, fb = component(za), component(zb)
            all(isfinite, (fa, fb)) && !iszero(fa) && !iszero(fb) &&
                signbit(fa) != signbit(fb) || continue
            t = ComplexRegions.brent(u -> component(point(side, u) - target),
                zero(T), one(T), T(4) * eps(T))
            isfinite(t) && zero(T) <= t <= one(T) || continue
            sr = _clean_param(d, T(k - 1) + T(t))
            matches(sr) && return sr
        end
    end

    throw(ArgumentError("Cannot locate the reflected point for curve parameter $s; the mesh cannot be refined at both nodes"))
end

# %% _mesh_params
# Return a sorted mesh retaining endpoints; specialize its arithmetic on T.
function _mesh_params(d::D, params, extra, ::Type{T}) where {D,T}
    L = T(length(d))
    tol = T(100eps(Float64)) * max(one(T), abs(L))

    # Start with both endpoints and collect cleaned interior parameters.
    vals = T[zero(T), L]
    for collection in (params, extra), s in collection
        value = _clean_param(d, T(s))
        (abs(value) <= tol || abs(value - L) <= tol) && continue
        push!(vals, value)
    end

    # Merge nearby parameters in sorted order.
    sort!(vals)
    out = T[]
    for value in vals
        clean = abs(value) <= tol ? zero(T) : abs(value - L) <= tol ? L : value
        if isempty(out) ||
                abs(clean - last(out)) > tol * max(one(T), abs(clean), abs(last(out)))
            push!(out, clean)
        end
    end

    # Keep both parameter endpoints even when they represent the same curve point.
    first(out) != zero(T) && pushfirst!(out, zero(T))
    last(out) != L && push!(out, L)
    return out
end

# %% _same
# Compare two nodes or values with a scaled tolerance.
_same(z, w) = abs(z - w) <= 100eps(Float64) * max(1, abs(z), abs(w))

# %% _node_hash_key
# Return the integer lookup key, specializing on the point and cell-width types.
_node_hash_key(z, h) = (round(Int, Float64(real(z)) / h), round(Int, Float64(imag(z)) / h))

# %% _lookup_node_index
# Return the exact or nearest matching node index, or nothing when none is found.
function _lookup_node_index(table, zgrid, z, h)
    k0 = _node_hash_key(z, h)
    for idx in get(table, k0, ())
        zgrid[idx] == z && return idx
    end
    nearest, distance = nothing, Inf

    # The distinctness tolerance spans up to two neighboring hash cells.
    for di in -2:2, dj in -2:2, idx in get(table, (k0[1] + di, k0[2] + dj), ())
        δ = abs(zgrid[idx] - z)
        if _same(zgrid[idx], z) && δ < distance
            nearest, distance = idx, δ
        end
    end
    return nearest
end

# %% _projective_reduce!
# Update homogeneous numerator/denominator arrays through one block.
# Common normalization preserves each ratio during construction.
function _projective_reduce!(hnum, hden, a, b, zgrid, pair)
    # Hoist the affine center; keep normalized construction separate from evaluation.
    affine = length(b) == 2
    midpoint = affine ? oftype(float(real(pair[1])), 0.5) * pair[1] +
        oftype(float(real(pair[1])), 0.5) * pair[2] : zero(eltype(b))
    center_value = affine ? b[1] + b[2] * midpoint : zero(eltype(b))

    # Apply h -> a/(h-b) to homogeneous pairs without dividing by h's denominator.
    @inbounds for k in eachindex(zgrid)
        N = hnum[k]
        D = hden[k]
        q = _evaluate_block(a, pair, zgrid[k]; numerator=true)
        bj = affine ? center_value + b[2] * (zgrid[k] - midpoint) : _evaluate_block(b, pair, zgrid[k])
        hnum[k] = q * D
        hden[k] = N - bj * D

        # Normalize numerator and denominator together to limit their growth.
        s = max(abs(hnum[k]), abs(hden[k]))
        if isfinite(s) && s > 0
            hnum[k] /= s
            hden[k] /= s
        end
    end
    return hnum, hden
end

# %% _reduced_values
# Return reduced sample values, reusing workspace when supplied.
function _reduced_values(ztest, ftest, bpolys, apolys, pairs; workspace=nothing)
    # Reuse compatible storage and initialize each value as a homogeneous pair.
    S = promote_type(ComplexF64, eltype(ztest), eltype(ftest))
    hnum, hden = workspace === nothing || eltype(workspace[1]) !== S ?
        (Vector{S}(undef, length(ftest)), Vector{S}(undef, length(ftest))) : workspace
    resize!(hnum, length(ftest))
    resize!(hden, length(ftest))
    hnum .= ftest
    fill!(hden, one(S))

    # Encode infinite input values without storing an infinite numerator.
    for k in eachindex(hnum)
        if isinf(hnum[k]) && !isnan(hnum[k])
            hnum[k], hden[k] = one(S), zero(S)
        end
    end

    # Reduce through all blocks before converting the final pairs back to values.
    for j in eachindex(bpolys)
        _projective_reduce!(hnum, hden, apolys[j], bpolys[j], ztest, pairs[j])
    end

    return hnum .= _qtcf_ratio.(hnum, hden)
end

# %% _reflected_continuum_samples
# Return unique sample points, values, and parameters including their reflections.
function _reflected_continuum_samples(f::F, d::D, zgrid, fgrid, sgrid, reflection::R) where {F,D,R}
    # Use one spatial lookup for original and reflected samples.
    S = promote_type(ComplexF64, eltype(zgrid), eltype(fgrid))
    T = promote_type(typeof(float(sgrid[firstindex(sgrid)])), typeof(float(length(d))))
    zout, fout, sout = S[], S[], T[]
    h = 64eps(Float64) * maximum(z -> Float64(abs(z)), S.(zgrid); init=1.0)
    table = Dict{Tuple{Int,Int},Vector{Int}}()

    for k in eachindex(zgrid), reflected in (false, true)
        # Reuse original values; locate and evaluate reflected points on the curve.
        s = T(sgrid[k])
        if reflected
            s = _reflection_param(d, s, reflection)
            (iszero(s) || s == length(d)) && continue
            z = S(point(d, s))
            y = S(f(z))
        else
            z, y = S(zgrid[k]), S(fgrid[k])
        end

        # Record each point once, keeping its value and cleaned parameter together.
        _lookup_node_index(table, zout, z, h) !== nothing && continue
        push!(zout, z)
        push!(fout, y)
        push!(sout, T(_clean_param(d, s)))
        push!(get!(Vector{Int}, table, _node_hash_key(z, h)), length(zout))
    end
    return zout, fout, sout
end

# %% _qtcf_converged
# Test a recorded error and allowed poles, caching the poles when needed.
function _qtcf_converged(record, threshold, allowed)
    isfinite(record.error) && record.error <= threshold || return false
    allowed === true && return true

    # Compute poles only after the error test passes and a pole check is required.
    ismissing(record.poles) && (record.poles = poles(record.interpolant))
    return all(allowed, record.poles)
end

# %% [markdown]
# ## 6. Continuum approximation, discrete approximation, and symmetry adapters
#
# Continuum `approximate` samples and refines a curve; discrete `approximate`
# uses a fixed grid. Their local block acceptors retain QTCF's finite-block,
# symmetry, and recording policies. Both reuse RFA histories, status types,
# and best-result selection. The two symmetry adapters unwrap the maps before
# calling these fitting methods. `check` includes the whole discrete grid.

# %% approximate — continuum
"""
    approximate(f, domain, QTCF(); tol, max_degree=100, max_iter=max_degree, kwargs...)

Construct a continuum greedy QTCF, returning the package's approximation and
history wrapper. `refinement=3` and `initial_refinement=15` control sampling;
`stagnation=5` limits unsuccessful block updates. Budgets count blocks, not
individual nodes. `allowed=true`, a pole predicate, or `allowed=:strict`
selects the same pole policy as TCF. See [`symmetric`](@ref) for pairing.
The mesh spans the full parameter interval, but test samples lie inside its
subintervals, as in TCF. The first block uses `point(domain, 0)` and its distinct
reflection, when supplied; these selected endpoints are not test samples.

New blocks require finite reduced values and coefficients. Every convergent is
recorded; error checks include selected nodes, which cannot be selected again.
Without a reflection, the first block
contains one node and subsequent blocks select the two largest residuals.
"""
function approximate(
    f::Function, d::ComplexCurveOrPath, ::QuadraticThiele;
    float_type::Type{F} = promote_type(real_type(d), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed = true,
    max_degree::Int = 100,
    max_iter::Int = max_degree,
    refinement::Int = 3,
    initial_refinement::Int = 15,
    stagnation::Int = 5,
    reflection::Union{Nothing,Function} = nothing,
    value_reflection::Union{Nothing,Function} = nothing
    ) where {F}

    # Resolve pole and symmetry options.
    if allowed == :strict
        # Only allow poles off the curve.
        allowed = z -> ComplexRegions.dist(z, d) > tol
    end
    reflection === nothing && value_reflection !== nothing && throw(ArgumentError(
        "value_reflection requires a reflection map"))

    # Initialize the mesh and reusable construction storage.
    paired = reflection !== nothing
    Tpar, S = typeof(float(length(d))), complex(float_type)
    L = Tpar(length(d))
    num_ref = max(refinement, initial_refinement)
    # Start with [0, 1] on a curve; multi-curve paths use their full native interval.
    s_nodes = Tpar[zero(Tpar), L]
    maxpoints = max(50, paired ? 2max_iter * num_ref + length(s_nodes) + 20 : 2max_iter + 20)
    path = DiscretizedPath(d, s_nodes; refinement=num_ref, maxpoints)
    sample_f = S[]

    # Return mesh points, sampled values, and parameters using a reusable buffer.
    function sample_path()
        indices = CartesianIndices((1:length(s_nodes)-1, 2:num_ref+1))
        z, s = vec(view(path.points, indices)), vec(view(path.params, indices))
        resize!(sample_f, length(z))
        _sweep!(sample_f, S ∘ f, z, eachindex(z))
        return z, sample_f, s
    end

    ztest, ftest, stest = paired ? sample_path() : (S[], S[], Tpar[])
    R = promote_type(ComplexF64, eltype(ztest), eltype(ftest))
    reduced_workspace = (R[], R[])
    fmax = paired ? maximum(abs, ftest) : zero(float_type)
    active_h = 64eps(Float64) * maximum(z -> Float64(abs(z)),
        paired ? ztest : S[point(d, zero(Tpar)), point(d, L)]; init=1.0)
    active_table = Dict{Tuple{Int,Int},Vector{Int}}()
    pairs = Tuple{S,S}[]
    bpolys, apolys, zact, yact = Vector{S}[], Vector{S}[], S[], S[]
    parity_fixed_index, parity_center = 0, zero(S)
    value_mode = :generic

    # Check value symmetry and detect equal values at reflected points.
    if paired
        probe_values, reflected_values = S[], S[]
        for sprobe in stest
            zprobe = S(point(d, _clean_param(d, sprobe)))
            zrprobe = S(reflection(zprobe))
            _same(zprobe, zrprobe) && value_reflection === nothing && continue
            yprobe, yrprobe = S(f(zprobe)), S(f(zrprobe))
            if value_reflection !== nothing
                mapped = S(value_reflection(zprobe, yprobe, zrprobe))
                isapprox(mapped, yrprobe; rtol=1e-10, atol=1e-12) || throw(ArgumentError(
                    "value_reflection does not match f at reflected continuum points"))
                yrprobe = mapped
            end
            _same(zprobe, zrprobe) && continue
            push!(probe_values, yprobe)
            push!(reflected_values, yrprobe)
            length(probe_values) == 16 && break
        end

        if !isempty(probe_values)
            V = promote_type(ComplexF64, eltype(probe_values), eltype(reflected_values))
            y, yr = V.(probe_values), V.(reflected_values)
            scale = max(1.0, maximum(abs, y; init=0.0), maximum(abs, yr; init=0.0))
            equality = maximum(abs, yr .- y; init=0.0) / scale
            value_mode = equality <= 1000eps(Float64) ? :equality : :generic
        end
    end

    # Try a block at the supplied parameters; return whether it was accepted.
    function accept_block!(s1raw, s2raw=s1raw)
        # Locate the proposed node and its optional partner.
        s1 = _clean_param(d, paired ? s1raw : Tpar(s1raw))
        z1 = S(point(d, _clean_param(d, s1)))
        s2 = paired ? _reflection_param(d, s1, reflection) :
            _clean_param(d, Tpar(s2raw))
        z2 = S(paired ? reflection(z1) : point(d, _clean_param(d, s2)))

        # Reject reuse of either node, including coincident curve endpoints.
        any(z -> _same(z, z1) || _same(z, z2), zact) && return false
        fixed = _same(z1, z2)

        # Evaluate and reduce the data before constructing the block.
        y1 = S(f(z1))
        y2 = paired && fixed ? y1 : paired && value_reflection !== nothing ?
            S(value_reflection(z1, y1, z2)) : S(f(z2))
        block_reduced = isempty(bpolys) ? (y1, y2) : paired && fixed ?
            _reduced_values((z1,), (y1,), bpolys, apolys, pairs; workspace=reduced_workspace) :
            _reduced_values((z1, z2), (y1, y2), bpolys, apolys, pairs; workspace=reduced_workspace)
        h1, h2 = first(block_reduced), last(block_reduced)
        all(isfinite, (y1, y2, h1, h2)) || return false
        if paired
            fmax = max(fmax, abs(y1), abs(y2))
        end

        # Fit the denominator, retaining odd reduced parity about the fixed point.
        b = convert(Vector{S}, _block_line(z1, h1, z2, h2; symmetric=paired))
        if paired && parity_fixed_index > 0 && !fixed &&
                _same((z1 + z2) / S(2), parity_center) &&
                Float64(abs(-h2 - h1)) <= 100eps(Float64) * max(1.0, Float64(abs(h1)), Float64(abs(-h2)))
            slope = (h2 - h1) / (z2 - z1)
            midpoint = (z1 + z2) / S(2)
            b = S[-slope * midpoint, slope]
        end
        all(isfinite, b) || return false
        a = convert(Vector{S}, _factor_poly(z1, z2))
        if paired && fixed && value_mode === :equality && parity_fixed_index == 0
            parity_fixed_index, parity_center = length(bpolys) + 1, z1
        end

        # Store the accepted block and mark its nodes as selected.
        push!(pairs, (z1, z2))
        push!(bpolys, b)
        push!(apolys, a)
        for (z, y) in (fixed ? ((z1, y1),) : ((z1, y1), (z2, y2)))
            push!(zact, z)
            push!(yact, y)
            push!(get!(Vector{Int}, active_table, _node_hash_key(z, active_h)), length(zact))
        end

        # Refine at both parameters.
        if paired
            s_nodes = _mesh_params(d, s_nodes, (s1, s2), typeof(s1))
        else
            fmax = max(fmax, abs(y1), abs(y2))
            s_nodes = _mesh_params(d, s_nodes, (s1, s2),
                promote_type(eltype(s_nodes), typeof(s1), typeof(s2)))
        end
        return true
    end

    # Initialize rational approximation and history.
    accept_block!(zero(Tpar)) || error("Initial QTCF continuum node was not admissible.")
    r = _qtcf_from_blocks(pairs, bpolys, apolys, zact, yact)
    history = [IterationRecord(r, NaN, missing)]
    best_err, stop_reason = Inf, :max_degree
    blocks = length(bpolys)
    last_improve = paired ? 0 : blocks
    reset!(path, s_nodes; refinement=num_ref)

    # Main iteration
    while true
        # Sample the mesh.
        ztest, ftest, stest = sample_path()
        fmax = max(fmax, maximum(abs, ftest))

        # Complete reflected samples and assess every selected node as well.
        if paired
            ztest, ftest, stest = _reflected_continuum_samples(f, d, ztest, ftest, stest, reflection)
        end
        err = abs.(ftest .- r.(ztest))
        node_error = maximum(zip(zact, yact); init=0.0) do (z, y)
            e = abs(y - r(z))
            isfinite(e) ? e : Inf
        end
        err_max = max(maximum(err; init=0.0), node_error)

        # Selected nodes contribute to the error, but cannot enter another block.
        for k in eachindex(ztest)
            _lookup_node_index(active_table, zact, paired ? ztest[k] : S(ztest[k]), active_h) !== nothing &&
                (err[k] = -Inf)
        end

        # Record every stage and apply the tolerance and block limits.
        history[end].error = err_max
        if err_max < best_err
            err_max < 0.999 * best_err && (last_improve = blocks)
            best_err = err_max
        end
        if _qtcf_converged(history[end], tol * fmax, allowed)
            stop_reason = :converged
            break
        end
        blocks >= max_iter && break
        recent = @view history[max(1, end-stagnation+1):end]
        if blocks - last_improve >= stagnation ||
                (length(history) >= stagnation + 5 && minimum(h.error for h in recent) > 5best_err)
            stop_reason = :stagnated
            break
        end

        # Attempt the highest-scoring node with its prescribed or greedy partner.
        accepted = false
        _, i = findmax(err)
        if err[i] != -Inf
            if paired
                accepted = accept_block!(stest[i])
            else
                second = copy(err)
                for k in eachindex(second)
                    _same(ztest[k], ztest[i]) && (second[k] = -Inf)
                end
                _, j = findmax(second)
                second[j] != -Inf && (accepted = accept_block!(stest[i], stest[j]))
            end
        end
        if !accepted
            stop_reason = :node_failure
            break
        end

        # Save the new stage and refine the mesh for the next iteration.
        r = _qtcf_from_blocks(pairs, bpolys, apolys, zact, yact)
        push!(history, IterationRecord(r, NaN, missing))
        blocks = length(bpolys)
        num_ref > refinement && (num_ref = max(refinement, num_ref - 1))
        try
            reset!(path, s_nodes; refinement=num_ref)
        catch
            stop_reason = :refinement
            break
        end
    end

    # Return the best acceptable approximation and its history.
    assessed = findall(record -> !isnan(record.error), history)
    best = isempty(assessed) ? lastindex(history) :
        assessed[best_acceptable(view(history, assessed), allowed)]
    best_r = copy(history[best].interpolant)
    stop = ConvergenceStatus(stop_reason, best, history)
    return ContinuumApproximation(f, d, best_r, allowed, path, history, stop)
end

# %% approximate — discrete
"""
    approximate(y, z, QTCF(); max_iter=length(y), tol, stagnation=5, kwargs...)

Fit sampled values `y` on the grid `z` with quadratic Thiele blocks.
`reflection` optionally pairs nodes on a reflection-closed grid; otherwise
the two largest residuals select each new pair. `max_iter` counts blocks. The stopping
error is relative to the sampled target scale.
Every convergent is recorded. Error checks use the whole grid, including selected
nodes; only unselected nodes are eligible for the next block.
For equal-value symmetry, alternating stages use the preceding paired error
to guide selection.
Use [`symmetric`](@ref) to declare point and function-value symmetry.
"""
function approximate(
    y::AbstractVector{T}, z::AbstractVector{Z}, ::QuadraticThiele;
    float_type::Type{F} = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = length(y),
    stagnation::Int = 5,
    reflection::Union{Nothing,Function} = nothing
    ) where {T<:Number,Z<:Number,F}

    # Choose the scalar TCF startup node.
    _, idx_min = findmin(abs, y)
    zgrid, fgrid = z, y
    tol = convert(float_type, tol)
    paired = reflection !== nothing

    # Build the node lookup and locate each reflected partner.
    S = complex(float_type)
    zgrid, fgrid = S.(zgrid), S.(fgrid)
    used = zeros(Int, length(zgrid))
    hash_step = 64eps(Float64) * maximum(z -> Float64(abs(z)), zgrid; init=1.0)
    table = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in eachindex(zgrid)
        push!(get!(Vector{Int}, table, _node_hash_key(zgrid[i], hash_step)), i)
    end
    partner = paired ? [something(_lookup_node_index(table, zgrid,
        reflection(z), hash_step), 0) for z in zgrid] : Int[]
    any(iszero, partner) && throw(ArgumentError(
        "Each reflected point must be present in the discrete grid."))

    # Initialize the stopping threshold, blocks, and reduced-data buffers.
    fmax = maximum(abs, fgrid)
    stop_tol = tol * fmax
    pairs = Tuple{S,S}[]
    bpolys, apolys, zact, yact = Vector{S}[], Vector{S}[], S[], S[]
    hnum, hden = copy(fgrid), ones(S, length(fgrid))
    reference_error = float_type[]
    value_mode = :generic
    if paired
        for i in eachindex(zgrid)
            _same(zgrid[i], reflection(zgrid[i])) && continue
            yr = fgrid[partner[i]]
            if !(abs(yr - fgrid[i]) <= 1000eps(float_type) * max(1, fmax))
                value_mode = :generic
                break
            end
            value_mode = :equality
        end
    end
    parity_fixed_index, parity_center = 0, zero(S)

    # Try a block at grid indices i and j.
    function accept_pair!(i, j)
        (iszero(used[i]) && iszero(used[j])) || return false

        # Read target and reduced values at both grid points.
        z1, y1, h1 = zgrid[i], fgrid[i], isempty(bpolys) ? fgrid[i] : _qtcf_ratio(hnum[i], hden[i])
        (isfinite(y1) && isfinite(h1)) || return false
        z2, y2 = zgrid[j], fgrid[j]
        h2 = isempty(bpolys) ? fgrid[j] : _qtcf_ratio(hnum[j], hden[j])
        (isfinite(y2) && isfinite(h2)) || return false

        # Fit a finite block and append its coefficients and interpolation data.
        b = convert(Vector{S}, _block_line(z1, h1, z2, h2; symmetric=paired))
        fixed = _same(z1, z2)
        if paired && parity_fixed_index > 0 && !fixed &&
                _same((z1 + z2) / S(2), parity_center) &&
                abs(h1 + h2) <= 100eps(float_type) * max(abs(h1), abs(h2))
            slope = (h2 - h1) / (z2 - z1)
            b = S[-slope * ((z1 + z2) / S(2)), slope]
        end
        all(isfinite, b) || return false
        a = convert(Vector{S}, _factor_poly(z1, z2))
        if paired && fixed && value_mode === :equality && parity_fixed_index == 0
            parity_fixed_index, parity_center = length(bpolys) + 1, z1
        end
        push!(pairs, (z1, z2))
        push!(bpolys, b)
        push!(apolys, a)
        for (z, y) in (fixed ? ((z1, y1),) : ((z1, y1), (z2, y2)))
            push!(zact, z)
            push!(yact, y)
        end

        # Exclude selected-node neighbors using the continuum distinctness rule.
        used .= ifelse.(iszero.(used) .& (_same.(zgrid, z1) .| _same.(zgrid, z2)), length(bpolys), used)
        _projective_reduce!(hnum, hden, a, b, zgrid, (z1, z2))
        return true
    end

    # Start at the point of minimum target magnitude.
    seed = idx_min
    accept_pair!(seed, paired ? partner[seed] : seed) ||
        error("Initial QTCF node matching the scalar TCF start was not admissible.")

    # Initialize rational approximation and history; retain node-selection marks.
    r = _qtcf_from_blocks(pairs, bpolys, apolys, zact, yact)
    history = IterationRecord{typeof(r),real_type(S),S}[]
    stop_reason, best_err = :max_degree, Inf
    last_improve = length(bpolys)

    # Main iteration
    while true
        # Measure the whole-grid error before excluding selected candidates.
        err = abs.(fgrid .- r.(zgrid))
        err_max = maximum(eachindex(err); init=0.0) do i
            !iszero(used[i]) && !isfinite(err[i]) ? Inf : err[i]
        end
        err .= ifelse.(iszero.(used), err, -Inf)
        blocks = length(bpolys)
        push!(history, IterationRecord(r, err_max, missing))
        paired_stage = value_mode !== :equality || parity_fixed_index == 0 ||
            blocks <= parity_fixed_index || iseven(blocks - parity_fixed_index)
        if paired_stage && value_mode === :equality && parity_fixed_index > 0
            resize!(reference_error, length(err))
            copyto!(reference_error, err)
        end
        if isfinite(err_max) && err_max < best_err
            err_max < 0.999 * best_err && (last_improve = blocks)
            best_err = err_max
        end

        # Stop on stagnation, convergence, or the block limit.
        if length(history) > 1 && (blocks - last_improve >= stagnation ||
                (length(history) >= stagnation + 5 &&
                    minimum(h.error for h in history[end-stagnation+1:end]) > 5best_err))
            stop_reason = :stagnated
            break
        end
        blocks < max_iter && !_qtcf_converged(history[end], stop_tol, allowed) || break

        # Combine current and preceding paired errors between paired stages.
        score = err
        if !paired_stage
            score = ifelse.(isfinite.(err) .& isfinite.(reference_error),
                sqrt.(max.(err, 0)) .* sqrt.(max.(reference_error, 0)), err)
            iszero(maximum(score; init=zero(float_type))) && (score = err)
        end

        # Attempt the highest-scoring block and update the approximant.
        accepted = false
        _, i = findmax(score)
        if score[i] != -Inf
            if paired
                accepted = accept_pair!(i, partner[i])
            else
                second = copy(score)
                for k in eachindex(second)
                    _same(zgrid[k], zgrid[i]) && (second[k] = -Inf)
                end
                _, j = findmax(second)
                second[j] != -Inf && (accepted = accept_pair!(i, j))
            end
        end
        if !accepted
            stop_reason = all(!iszero, used) ? :exhausted : :node_failure
            break
        end
        r = _qtcf_from_blocks(pairs, bpolys, apolys, zact, yact)
    end

    # Return the best acceptable approximation and its remaining-data mask.
    _qtcf_converged(history[end], stop_tol, allowed) && (stop_reason = :converged)
    assessed = findall(record -> !isnan(record.error), history)
    best = isempty(assessed) ? lastindex(history) :
        assessed[best_acceptable(view(history, assessed), allowed)]
    best_r = copy(history[best].interpolant)
    stop = ConvergenceStatus(stop_reason, best, history)
    return DiscreteApproximation(y, z, best_r, iszero.(used) .| (used .> best), allowed, history, stop)
end

# %% check — complete discrete grid
"Check the QTCF error at every supplied grid point, including selected nodes."
function check(r::DiscreteApproximation{T,S,R}; quiet=false, verbose=!quiet) where {T,S,R<:QuadraticThiele}
    z = r.domain
    err = r.data - r.(z)
    verbose && @info f"Max error is {norm(err, Inf):.2e}"
    return z, err
end

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
        # Validate sampled partner values before dispatching to the discrete method.
        y = f.domain.(d)
        for (z, value) in zip(d, y)
            isapprox(f.reflection(value), f.domain(reflection(z));
                rtol=1e-10, atol=1e-12) || throw(ArgumentError(
                "The function reflection does not match f at reflected sample points"))
        end
        return approximate(y, d, QuadraticThiele(); reflection, kwargs...)
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
