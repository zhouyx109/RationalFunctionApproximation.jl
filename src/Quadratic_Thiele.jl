# -*- coding: utf-8 -*-
# %% [markdown]
# # Quadratic Thiele continued fractions
#
# ## 1. Types and accessors
#
# `_QTCFBlocks` stores blocks; `QuadraticThiele` uses Thiele storage for singletons.
# `_QTCFChart` retains normalized polynomial factors at a selected common zero.
# `nodes`, `values`, `degrees`, and `degree` support RFA results and display.
# `convert` changes numeric type; `copy` saves snapshots; `isreal` checks coefficients.
# `getproperty` builds continuants for polynomial queries.

# %%
import Base: values

export QuadraticThiele, QTCF, Symmetric, symmetric

# %% _QTCFGroup — normalized polynomial transfer coefficients
mutable struct _QTCFGroup{S}
    head::Vector{S}
    gamma::Vector{S}
    center::S
    complement::NTuple{3,S}
    ratio::S
    scale::S
end

struct _QTCFChart{S}
    groups::Vector{_QTCFGroup{S}}
end

Base.copy(g::_QTCFGroup{S}) where {S} = _QTCFGroup{S}(
    copy(g.head), copy(g.gamma), g.center, g.complement, g.ratio, g.scale)
Base.copy(c::_QTCFChart{S}) where {S} = _QTCFChart{S}(copy.(c.groups))

_convert_chart(::Type{S}, chart) where {S} = chart === nothing ? nothing :
    _QTCFChart{S}([_QTCFGroup{S}(S.(g.head), S.(g.gamma), S(g.center),
        S.(g.complement), S(g.ratio), S(g.scale)) for g in chart.groups])

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
vectors, one per block. Singleton storage delegates to Thiele. When the first
pair introduces a selected common zero, computed blocks use normalized polynomial
transfer factors; expanded denominator weights remain available. Continuants are lazy.
"""
mutable struct QuadraticThiele{T,S} <: AbstractRationalInterpolant{T,S}
    data::Union{Thiele{T,S},_QTCFBlocks{S}}
    polynomials::Union{Nothing,Tuple{Vector{S},Vector{S}}}
    chart::Union{Nothing,_QTCFChart{S}}

    function QuadraticThiele{T,S}(data, polynomials, chart=nothing) where {T<:AbstractFloat,S<:RealComplex{T}}
        return new{T,S}(data, polynomials, chart)
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
    chart = _convert_chart(V, r.chart)
    return QuadraticThiele{F,V}(data, nothing, chart)
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
    return QuadraticThiele{T,S}(data, p === nothing ? nothing : copy.(p),
        r.chart === nothing ? nothing : copy(r.chart))
end

# %% Base.getproperty
# Return stored nodes/values or build polynomial continuants.
@inline function Base.getproperty(r::QuadraticThiele{T,S}, name::Symbol) where {T,S}
    d = getfield(r, :data)
    if name === :P || name === :Q
        p = getfield(r, :polynomials)
        if p === nothing && r.chart !== nothing
            r.polynomials = p = _chart_polynomials(r)
        end
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
    name in (:data, :polynomials, :chart) && return getfield(r, name)
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
# `_chart_numden`, `_chart_evaluate!`, and `_chart_derivative!` evaluate the same
# normalized transfer factors; `_chart_local_pair` cancels a selected linear factor.
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
    d isa Thiele && !isfinite(z) &&
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

    S = promote_type(eltype(r), typeof(z))
    value = d isa Thiele ?
        convert(S, evaluate(d, z, method)) : _evaluate(r, z, method)

    # Cancel a selected common factor before evaluating the quotient.
    if !isfinite(value) && any(==(z), nodes(r))
        cancelled = _evaluate_common_zero(r, z)
        cancelled === nothing || return cancelled
    end
    return value
end

# %% evaluate — array
# Promote stored coefficients with the query type.
function evaluate(r::QuadraticThiele, z::AbstractArray{<:Number})
    t = similar(z, promote_type(eltype(r), eltype(z)))
    return evaluate!(t, r, z)
end

# %% evaluate! — allocate scratch
"""
    evaluate!(t, r::QuadraticThiele, z, [a, b])

Write evaluations into `t` using matching, nonoverlapping scratch arrays.
The OneDiv array recurrence omits scalar underflow rescaling, as in Thiele.
Selected common zeros use the same cancellation as scalar evaluation.
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

    # Use the array pair recurrence without scalar underflow checks.
    if r.chart !== nothing
        _chart_evaluate!(t, r, z, a, b)
    elseif d isa Thiele
        evaluate!(t, d, z, a, b; fallback)
    else
        m = length(d.bpolys)
        if m == 1
            p, pair = d.bpolys[1], d.pairs[1]
            if length(p) == 1
                t .= p[1]
            else
                half = oftype(float(real(pair[1])), 0.5)
                midpoint = half * pair[1] + half * pair[2]
                value = p[1] + p[2] * midpoint
                t .= value .+ p[2] .* (z .- midpoint)
            end
        else
            p, pair = d.bpolys[m], d.pairs[m]
            if length(p) == 1
                a .= p[1]
            else
                half = oftype(float(real(pair[1])), 0.5)
                midpoint = half * pair[1] + half * pair[2]
                value = p[1] + p[2] * midpoint
                a .= value .+ p[2] .* (z .- midpoint)
            end
            p, pair = d.apolys[m-1], d.pairs[m-1]
            if length(p) == 2
                b .= p[2] .* (z .- pair[1])
            else
                b .= (p[3] .* (z .- pair[1])) .* (z .- pair[2])
            end
            @inbounds for k in m-1:-1:2
                p, pair = d.bpolys[k], d.pairs[k]
                if length(p) == 1
                    b .= p[1] .* a .+ b
                else
                    half = oftype(float(real(pair[1])), 0.5)
                    midpoint = half * pair[1] + half * pair[2]
                    value = p[1] + p[2] * midpoint
                    b .= (value .+ p[2] .* (z .- midpoint)) .* a .+ b
                end
                p, pair = d.apolys[k-1], d.pairs[k-1]
                if length(p) == 2
                    a .*= p[2] .* (z .- pair[1])
                else
                    a .*= (p[3] .* (z .- pair[1])) .* (z .- pair[2])
                end
                a, b = b, a
            end
            p, pair = d.bpolys[1], d.pairs[1]
            if length(p) == 1
                t .= p[1] .+ b ./ a
            else
                half = oftype(float(real(pair[1])), 0.5)
                midpoint = half * pair[1] + half * pair[2]
                value = p[1] + p[2] * midpoint
                t .= (value .+ p[2] .* (z .- midpoint)) .+ b ./ a
            end
        end
    end

    # Resolve exceptional inputs and selected common zeros.
    @inbounds for i in eachindex(t, z)
        (isinf(z[i]) || (!isfinite(t[i]) && any(==(z[i]), nodes(r))) ||
            (fallback && (isnan(z[i]) || !isfinite(t[i])))) &&
            (t[i] = evaluate(r, z[i], OneDiv()))
    end
    return t
end

# %% _evaluate_classic
# Evaluate by the backward continued-fraction recurrence.
function _evaluate_classic(r::QuadraticThiele, z::Number)
    r.chart === nothing || return /( _chart_numden(r, z; classic=true)...)
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
    r.chart === nothing || return _chart_numden(r, z; underflow_guard)
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

# Evaluate the normalized head, tail coefficient, and anchored complement.
@inline _chart_head(g, z) = length(g.head) == 1 ? g.head[1] : g.head[1] + g.head[2]*z
@inline _chart_gamma(g, z) = length(g.gamma) == 1 ? g.gamma[1] : g.gamma[1] + g.gamma[2]*z
@inline _chart_complement(g, z) = evalpoly(z - g.center, g.complement)

# Return the final local pair and apply preceding polynomial transfer matrices.
function _chart_numden(r::R, z; underflow_guard=true, classic=false) where {R<:QuadraticThiele}
    groups, pairs = r.chart.groups, r.pairs
    g = last(groups)
    complete = !isempty(g.gamma)
    P = complete ? _chart_complement(g, z) : _chart_head(g, z)
    Q = complete ? _chart_gamma(g, z) : one(P)
    classic && ((P, Q) = (_chart_divide(P, Q), one(Q)))
    for k in length(groups)-1:-1:1
        g = groups[k]
        pair = pairs[2k]
        a = pair[1] == pair[2] ? z - pair[1] : (z - pair[1]) * (z - pair[2])
        if classic && isinf(P)
            P, Q = _chart_divide(_chart_complement(g, z), _chart_gamma(g, z)), one(Q)
            continue
        end
        P, Q = g.ratio * _chart_complement(g, z) * P + _chart_head(g, z) * a * Q,
            g.ratio * _chart_gamma(g, z) * P + a * Q
        if classic
            P, Q = P / Q, one(Q)
        elseif k > 1 && underflow_guard && guard(P)
            P *= 1e100
            Q *= 1e100
        end
    end
    return P, Q
end

# Apply the OneDiv transfer chain with the native unscaled array policy.
function _chart_evaluate!(out, r::R, z, a, b) where {R<:QuadraticThiele}
    groups, pairs = r.chart.groups, r.pairs
    g = last(groups)
    if isempty(g.gamma)
        h0, h1 = g.head[1], length(g.head)>1 ? g.head[2] : zero(eltype(a))
        a .= h0 .+ h1 .* z
        b .= one(eltype(b))
    else
        c0, c1, c2 = g.complement
        gamma0, gamma1 = g.gamma[1], length(g.gamma)>1 ? g.gamma[2] : zero(eltype(a))
        a .= c0 .+ (z .- g.center) .* (c1 .+ c2 .* (z .- g.center))
        b .= gamma0 .+ gamma1 .* z
    end
    temp = out
    for k in length(groups)-1:-1:1
        g, pair = groups[k], pairs[2k]
        c0, c1, c2 = g.complement
        h0, h1 = g.head[1], length(g.head)>1 ? g.head[2] : zero(eltype(a))
        gamma0, gamma1 = g.gamma[1], length(g.gamma)>1 ? g.gamma[2] : zero(eltype(a))
        if pair[1] == pair[2]
            temp .= g.ratio .* (c0 .+ (z .- g.center) .* (c1 .+ c2 .* (z .- g.center))) .* a .+
                (h0 .+ h1 .* z) .* (z .- pair[1]) .* b
            b .= g.ratio .* (gamma0 .+ gamma1 .* z) .* a .+ (z .- pair[1]) .* b
        else
            temp .= g.ratio .* (c0 .+ (z .- g.center) .* (c1 .+ c2 .* (z .- g.center))) .* a .+
                (h0 .+ h1 .* z) .* ((z .- pair[1]) .* (z .- pair[2])) .* b
            b .= g.ratio .* (gamma0 .+ gamma1 .* z) .* a .+
                ((z .- pair[1]) .* (z .- pair[2])) .* b
        end
        a, temp = temp, a
    end
    out .= a ./ b
    return out
end

# Differentiate the same polynomial transfer chain into caller-owned workspaces.
function _chart_derivative!(z, r::R, order, A, B;
    underflow_guard=true, normalize=false) where {R<:QuadraticThiele}
    groups, pairs = r.chart.groups, r.pairs
    g = last(groups)
    complete = !isempty(g.gamma)
    fill!(A, 0); fill!(B, 0)
    if complete
        A[1], B[1] = _chart_complement(g, z), _chart_gamma(g, z)
        order >= 1 && (A[2] = g.complement[2] + 2g.complement[3] * (z - g.center);
            B[2] = length(g.gamma) > 1 ? g.gamma[2] : zero(z))
        order >= 2 && (A[3] = 2g.complement[3])
    else
        A[1], B[1] = _chart_head(g, z), one(z)
        order >= 1 && length(g.head) > 1 && (A[2] = g.head[2])
    end
    if normalize
        scale = max(one(abs(z)), maximum(abs, A), maximum(abs, B))
        A ./= scale; B ./= scale
    end
    for k in length(groups)-1:-1:1
        g, pair = groups[k], pairs[2k]
        singleton = pair[1] == pair[2]
        a = singleton ? z - pair[1] : (z - pair[1]) * (z - pair[2])
        a1 = singleton ? one(z) : (z - pair[1]) + (z - pair[2])
        a2 = singleton ? zero(z) : oftype(z, 2)
        h, h1 = _chart_head(g, z), length(g.head) > 1 ? g.head[2] : zero(z)
        c, c1, c2 = _chart_complement(g, z),
            g.complement[2] + 2g.complement[3] * (z - g.center), 2g.complement[3]
        gamma, gamma1 = _chart_gamma(g, z), length(g.gamma) > 1 ? g.gamma[2] : zero(z)
        M11 = (g.ratio*c, g.ratio*c1, g.ratio*c2, zero(z))
        M12 = (h*a, h1*a+h*a1, 2h1*a1+h*a2, 3h1*a2)
        M21 = (g.ratio*gamma, g.ratio*gamma1, zero(z), zero(z))
        M22 = (a, a1, a2, zero(z))
        for q in order:-1:0
            p, d = zero(eltype(A)), zero(eltype(B))
            for s in 0:min(q, 3)
                choose = binomial(q, s)
                p += choose * (M11[s+1] * A[q-s+1] + M12[s+1] * B[q-s+1])
                d += choose * (M21[s+1] * A[q-s+1] + M22[s+1] * B[q-s+1])
            end
            A[q+1], B[q+1] = p, d
        end
        if normalize
            scale = max(one(abs(z)), maximum(abs, A), maximum(abs, B))
            A ./= scale; B ./= scale
        end
        if k > 1 && underflow_guard && any(guard, A)
            A .*= 1e50
            B .*= 1e50
        end
    end
    return A, B
end

# Expand the anchored complement only for polynomial and pencil construction.
function _chart_coefficients(g)
    c0, c1, c2 = g.complement
    return [c0 - g.center*c1 + g.center^2*c2, c1 - 2g.center*c2, c2]
end

# Expand the transfer chain into continuants without cancelling common factors.
function _chart_polynomials(r::QuadraticThiele)
    groups, a = r.chart.groups, r.apolys
    S = eltype(r)
    g = last(groups)
    P, Q = isempty(g.gamma) ? (copy(g.head), S[1]) :
        (_chart_coefficients(g), copy(g.gamma))
    for k in length(groups)-1:-1:1
        g, ak = groups[k], a[2k]
        P, Q = _padd(_pmul(g.ratio .* _chart_coefficients(g), P),
                _pmul(_pmul(g.head, ak), Q)),
            _padd(_pmul(g.ratio .* g.gamma, P), _pmul(ak, Q))
    end
    return P, Q
end

# Retain constant and linear local coefficients for exact selected-factor cancellation.
function _chart_local_pair(r::R, z) where {R<:QuadraticThiele}
    groups, pairs = r.chart.groups, r.pairs
    g = last(groups)
    complete = !isempty(g.gamma)
    p0 = complete ? _chart_complement(g, z) : _chart_head(g, z)
    p1 = complete ? g.complement[2]+2g.complement[3]*(z-g.center) :
        length(g.head) > 1 ? g.head[2] : zero(z)
    q0 = complete ? _chart_gamma(g, z) : one(z)
    q1 = complete && length(g.gamma) > 1 ? g.gamma[2] : zero(z)
    for k in length(groups)-1:-1:1
        g, pair = groups[k], pairs[2k]
        a = _chart_a(pair, z)
        a1 = pair[1] == pair[2] ? one(z) : (z-pair[1])+(z-pair[2])
        h, h1 = _chart_head(g, z), length(g.head)>1 ? g.head[2] : zero(z)
        c, c1 = _chart_complement(g, z), g.complement[2]+2g.complement[3]*(z-g.center)
        gamma, gamma1 = _chart_gamma(g, z), length(g.gamma)>1 ? g.gamma[2] : zero(z)
        p0, p1, q0, q1 = g.ratio*c*p0+h*a*q0,
            g.ratio*(c*p1+c1*p0)+h*a*q1+(h1*a+h*a1)*q0,
            g.ratio*gamma*p0+a*q0,
            g.ratio*(gamma*p1+gamma1*p0)+a*q1+a1*q0
        if k > 1 && (guard(p0) || guard(p1))
            p0 *= 1e50; p1 *= 1e50; q0 *= 1e50; q1 *= 1e50
        end
    end
    return p0, p1, q0, q1
end

# %% _evaluate_common_zero
# Divide a selected common factor out of the local lower-tail polynomials.
function _evaluate_common_zero(r::QuadraticThiele, z)
    if r.chart !== nothing
        S = promote_type(eltype(r), typeof(z))
        p0, p1, q0, q1 = _chart_local_pair(r, S(z))
        iszero(p0) && iszero(q0) && all(isfinite, (p1, q1)) &&
            !(iszero(p1) && iszero(q1)) || return nothing
        return p1 / q1
    end
    pairs, b, a = r.pairs, r.bpolys, r.apolys
    j = findfirst(pair -> z == first(pair) || z == last(pair), pairs)
    m = length(b)
    (j === nothing || j == m) && return nothing
    S = promote_type(eltype(r), typeof(z))

    # Retain the constant and linear coefficients in u=t-z of the lower tail.
    p0 = S(_evaluate_block(b[m], pairs[m], z))
    p1 = length(b[m]) == 1 ? zero(S) : S(b[m][2])
    q0, q1 = one(S), zero(S)
    for k in m-1:-1:j+1
        b0 = _evaluate_block(b[k], pairs[k], z)
        b1 = length(b[k]) == 1 ? zero(S) : b[k][2]
        a0 = _evaluate_block(a[k], pairs[k], z; numerator=true)
        a1 = length(a[k]) == 2 ? a[k][2] :
            a[k][3] * ((z - pairs[k][1]) + (z - pairs[k][2]))
        p0, p1, q0, q1 = b0 * p0 + a0 * q0,
            b0 * p1 + b1 * p0 + (a0 * q1 + a1 * q0), p0, p1
        if guard(p0) || guard(p1)
            p0 *= 1e50
            p1 *= 1e50
            q0 *= 1e50
            q1 *= 1e50
        end
    end
    iszero(p0) && isfinite(q0) && !iszero(q0) && isfinite(p1) || return nothing

    # Divide the vanishing lower numerator and a_j by u, then resume backward.
    pair = pairs[j]
    factor = length(a[j]) == 2 ? a[j][2] :
        a[j][3] * (z - (z == first(pair) ? last(pair) : first(pair)))
    P, Q = _evaluate_block(b[j], pair, z) * p1 + factor * q0, p1
    for k in j-1:-1:1
        P, Q = _evaluate_block(b[k], pairs[k], z) * P +
            _evaluate_block(a[k], pairs[k], z; numerator=true) * Q, P
        if k > 1 && guard(P)
            P *= 1e100
            Q *= 1e100
        end
    end
    return all(isfinite, (P, Q)) && !(iszero(P) && iszero(Q)) ? P / Q : nothing
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

# %% _derivative! — workspace and quotient recovery
# Fill numerator/denominator derivatives and recover vals with one reciprocal.
function _derivative!(
    z, r::QuadraticThiele, order::Integer,
    A::AbstractVector, B::AbstractVector, vals::AbstractVector;
    underflow_guard=true, quotient=true, normalize=false
    )
    n = order + 1
    d, S = getfield(r, :data), promote_type(eltype(r), typeof(z))
    if r.chart !== nothing
        _chart_derivative!(z, r, order, A, B; underflow_guard, normalize)
    elseif d isa Thiele
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
        if normalize
            scale = max(one(abs(z)), maximum(abs, view(A, 1:n)),
                maximum(abs, view(B, 1:n)))
            view(A, 1:n) ./= scale; view(B, 1:n) ./= scale
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
            if normalize
                scale = max(one(abs(z)), maximum(abs, view(A, 1:n)),
                    maximum(abs, view(B, 1:n)))
                view(A, 1:n) ./= scale; view(B, 1:n) ./= scale
            end
            if k > 1 && underflow_guard && any(guard, view(A, 1:n))
                view(A, 1:n) .*= 1e50
                view(B, 1:n) .*= 1e50
            end
        end
    end
    quotient || return A, B

    # Divide the local Taylor coefficients by a selected simple common factor.
    if iszero(A[1]) && iszero(B[1]) && any(==(z), nodes(r))
        _derivative!(z, r, order+1, A, B, vals; underflow_guard, quotient=false)
        if all(isfinite, (A[2], B[2])) && !(iszero(A[2]) && iszero(B[2]))
            for q in 0:order
                A[q+1], B[q+1] = A[q+2]/(q+1), B[q+2]/(q+1)
            end
        end
    end

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
# Return commonly normalized numerator/denominator jets for residues.
function _evaluate_numden_derivs(r::QuadraticThiele, z::Number, space=nothing)
    d = getfield(r, :data)
    d isa Thiele && return _evaluate_numden_derivs(d, z)
    S = promote_type(eltype(r), typeof(z))
    A, B = space === nothing ?
        ntuple(_ -> Vector{S}(undef, 2), 2) : space
    _derivative!(z, r, 1, A, B, A; underflow_guard=false, quotient=false, normalize=true)
    return S(A[1]), S(B[1]), S(A[2]), S(B[2])
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
function derivative(r::QuadraticThiele, orders::AbstractVector{<:Integer})
    isempty(orders) && throw(ArgumentError("at least one derivative order is required"))
    all(>=(0), orders) || throw(ArgumentError("derivative orders must be nonnegative"))
    order, index = maximum(orders), 1 .+ orders
    A = Vector{complex(eltype(r))}(undef, order + 2)
    B = similar(A)
    vals = similar(A, order + 1)
    return function(z)
        _derivative!(z, r, order, A, B, vals)
        result = vals[index]
        return isreal(z) && isreal(r) ? real(result) : result
    end
end

# %% [markdown]
# ## 3. Poles, roots, residues, and partial fractions
#
# `roots` and `poles` assemble the numerator and denominator matrix pencils.
# Normalized transfer factors use their polynomial companion pencil.
# `cancel_common` optionally removes confirmed common factors from the pencil roots.
# `residues` reuses RFA's residue loop with `_evaluate_numden_derivs`.
# `convert(PartialFractions, r)` uses the shared root/residue matching construction.

# %% poles
"""
    poles(r::QuadraticThiele; cancel_common=false)

Return finite roots of the unreduced denominator pencil, as for Thiele.
Set `cancel_common=true` to remove confirmed selected common factors.
"""
poles(r::QuadraticThiele; cancel_common=false) = roots(r, 2; cancel_common)

# Build a regular companion pencil for a polynomial in the transfer representation.
function _chart_polynomial_pencil(p)
    n = something(findlast(!iszero, p), 1) - 1
    S = eltype(p)
    C, D = zeros(S, n, n), zeros(S, n, n)
    for j in 1:n-1
        C[j,j+1], D[j,j] = -one(S), one(S)
    end
    if n > 0
        C[n,:] .= view(p, 1:n)
        D[n,n] = p[n+1]
    end
    return C, D
end

# %% roots
"""
    roots(r::QuadraticThiele, first_block=1; cancel_common=false)

Return finite numerator (`first_block=1`) or denominator (`first_block=2`) roots.
Both cases use an unreduced matrix pencil. Normalized transfer factors use a
companion pencil of their continuants. With `cancel_common=true`, remove one root for each confirmed selected common factor.
"""
function roots(r::QuadraticThiele, first_block::Integer=1; cancel_common=false)
    data = getfield(r, :data)
    if r.chart !== nothing
        p = first_block == 1 ? r.P : r.Q
        C, D = _chart_polynomial_pencil(p)
        isempty(C) && return eltype(r.nodes)[]
        z = try
            filter(isfinite, eigvals(-C, D))
        catch
            _, _, _, _, alpha, beta = schur(complex(-C), complex(D))
            filter(isfinite, alpha ./ beta)
        end
        return cancel_common ? _cancel_common_roots!(z, r) : z
    elseif data isa Thiele
        z = first_block == 2 ? poles(data) : roots(data)
        return cancel_common ? _cancel_common_roots!(z, r) : z
    end
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

    z = try
        filter(isfinite, eigvals(-C, D))
    catch
        _, _, _, _, α, β = schur(complex(-C), complex(D))
        filter(isfinite, α ./ β)
    end
    return cancel_common ? _cancel_common_roots!(z, r) : z
end

# %% _cancel_common_roots!
# Remove the pencil eigenvalue associated with each confirmed common factor.
function _cancel_common_roots!(z, r::QuadraticThiele)
    for node in nodes(r)
        isempty(z) && break
        _evaluate_common_zero(r, node) === nothing && continue
        deleteat!(z, argmin(abs.(z .- node)))
    end
    return z
end

# %% residues
"""
    residues(r::QuadraticThiele; cancel_common=false, check_finite=false, radius=nothing)

Return poles and residues using the shared TCF residue loop.
"""
function residues(r::QuadraticThiele; cancel_common=false, check_finite=false, radius=nothing)
    z = poles(r; cancel_common)
    if isempty(z) || getfield(r, :data) isa Thiele
        return residues(r, z; check_finite, radius)
    end
    S = promote_type(eltype(r), eltype(z))
    space = (Vector{S}(undef, 2), Vector{S}(undef, 2))
    numden = (r, z) -> _evaluate_numden_derivs(r, z, space)
    return residues(r, z; check_finite, radius, numden)
end

# %% Base.convert — partial fractions
"""
    convert(PartialFractions, r::QuadraticThiele)

Convert a QTCF with simple computed poles to polynomial-plus-partial-fraction form.
"""
Base.convert(::Type{PartialFractions}, r::QuadraticThiele) = _partial_fractions(r)

_polynomial_degree(r::QuadraticThiele) = something(findlast(!iszero, r.P), 1) - 1

# %% [markdown]
# ## 4. Data construction, node updates, and blocks
#
# Constructors and `add_node!` reuse `_new_weight` for reduced values and
# `_qtcf_append_block!` for storage.
# The shared approximation engine calls `add_node!` for every sampling adapter.
# `_block_line` fits denominator coefficients; `_factor_poly` forms numerators.
# `_qtcf_from_blocks` also supports numeric conversion and arithmetic.

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
# Wrap owned arrays in scalar storage when possible, otherwise in block storage.
function _qtcf_from_blocks(pairs, bpolys, apolys, z::Vector{S}, y, chart=nothing) where {S}
    scalar = !isempty(z) && length(z) == length(pairs) == length(bpolys) == length(apolys) &&
        all(eachindex(z)) do j
            pairs[j][1] == pairs[j][2] == z[j] && length(bpolys[j]) == 1 &&
                isfinite(only(bpolys[j])) && length(apolys[j]) == 2 &&
                apolys[j][1] == -z[j] && apolys[j][2] == 1
        end

    data = scalar ? Thiele(z, y, S[only(b) for b in bpolys]) :
        _QTCFBlocks{S}(pairs, bpolys, apolys, z, y)
    return QuadraticThiele{real_type(S),S}(data, nothing, chart)
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
    chart = r.chart === nothing ? nothing : copy(r.chart)
    if chart !== nothing
        update = _chart_append!(chart, r.pairs, (z, z), (y, y), method)
        update === nothing && throw(ArgumentError("QTCF block weights must be finite"))
        h = first(first(update))
    end
    _qtcf_append_block!(r, (z,), (y,), (h,))
    r.chart = chart
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

    (_same(first(z), last(z)) || any(in(nodes(r)), z)) &&
        throw(ArgumentError("interpolation nodes must be distinct"))

    # Reduce both values against the current approximant before fitting the block.
    h = [_new_weight(method, r, z[k], y[k]) for k in eachindex(z)]
    any(isnan, h) && throw(NaNException("Adding QTCF block caused NaN reduced value"))
    all(isfinite, h) || throw(ArgumentError("QTCF block reduced values must be finite"))

    b = _block_line(first(z), first(h), last(z), last(h))
    chart = r.chart === nothing && length(r.bpolys) == 1 &&
        first(r.pairs)[1] == first(r.pairs)[2] &&
        iszero(_evaluate_block(b, (first(z), last(z)), first(r.nodes))) ? _chart_start(eltype(r), first(r.vals)) :
        r.chart === nothing ? nothing : copy(r.chart)
    if chart !== nothing
        update = _chart_append!(chart, r.pairs, (first(z), last(z)), (first(y), last(y)), method)
        update === nothing && throw(ArgumentError("QTCF block weights must be finite"))
        b = first(update)
    end
    _qtcf_append_block!(r, z, y, b)
    r.chart = chart
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

# Initialize the singleton head without changing its target value.
function _chart_start(::Type{S}, value) where {S}
    g = _QTCFGroup{S}(S[value], S[], zero(S),
        (zero(S), zero(S), zero(S)), one(S), one(S))
    return _QTCFChart{S}([g])
end

# Evaluate the monic numerator associated with a stored node block.
@inline _chart_a(pair, z) = pair[1] == pair[2] ?
    z - pair[1] : (z - pair[1]) * (z - pair[2])

# Use the native extended-complex convention at an exact zero denominator.
@inline function _chart_divide(n, d)
    S = promote_type(typeof(n), typeof(d))
    iszero(d) && return convert(S, iszero(n) ? NaN : Inf)
    isinf(n) && isfinite(d) && return convert(S, Inf)
    return n / d
end

# Use ordinary Classic reductions for the first two blocks.
function _chart_first_two(chart, pairs, z, y)
    first = chart.groups[1]
    u = _chart_divide(_chart_a(pairs[1], z), y - _chart_head(first, z))
    return _chart_divide(_chart_a(pairs[2], z), u - _chart_gamma(first, z))
end

# Return h_{2k-2}/S_{2k-1}, the input to a normalized odd head.
function _chart_input(chart, pairs, z, y, k)
    k == 1 && return y
    u = _chart_divide(_chart_first_two(chart, pairs, z, y), chart.groups[1].ratio)
    for index in 2:k-1
        g = chart.groups[index]
        a = _chart_a(pairs[2index], z)
        B = _chart_head(g, z)
        u = isfinite(u) ? _chart_divide(a * (u - B),
            g.ratio * (_chart_complement(g, z) - _chart_gamma(g, z) * u)) :
            _chart_divide(_chart_divide(a,
                _chart_divide(_chart_a(pairs[2index-1], z), u-B) - _chart_gamma(g, z)),
                g.ratio)
    end
    return u
end

# Apply prior groups to a homogeneous reduced-value pair.
function _chart_input_pair(chart, pairs, z, y, k)
    U, V = y, one(y)
    for index in 1:k-1
        g = chart.groups[index]
        a = _chart_a(pairs[2index], z)
        U, V = a * (U - _chart_head(g, z) * V),
            g.ratio * (_chart_complement(g, z) * V - _chart_gamma(g, z) * U)
    end
    return U, V
end

# Return the current reduction scaled by the last odd-head scale.
function _chart_last_reduced(chart, pairs, z, y, method::ThieleMethod=Classic())
    k = length(chart.groups)
    if method isa OneDiv
        g = chart.groups[k]
        U, V = _chart_input_pair(chart, pairs, z, y, k)
        return _chart_divide(_chart_a(pairs[2k], z) * (U - _chart_head(g, z) * V),
            _chart_complement(g, z) * V - _chart_gamma(g, z) * U)
    end
    k == 1 && return _chart_first_two(chart, pairs, z, y)
    g = chart.groups[k]
    u = _chart_input(chart, pairs, z, y, k)
    a = _chart_a(pairs[2k], z)
    B = _chart_head(g, z)
    return isfinite(u) ? _chart_divide(a * (u-B),
        _chart_complement(g, z) - _chart_gamma(g, z)*u) :
        _chart_divide(a, _chart_divide(_chart_a(pairs[2k-1], z), u-B) -
            _chart_gamma(g, z))
end

# Complete or append one group and return canonical coefficients and reduced data.
# Existing pairs do not yet include newpair. A rejected update leaves the chart intact.
function _chart_append!(chart::_QTCFChart{S}, pairs, newpair, targets,
    method::ThieleMethod=default_weight_method()) where {S}
    x, y = S(newpair[1]), S(newpair[2])
    fx, fy = S(targets[1]), S(targets[2])
    all(isfinite, (fx, fy)) || return nothing
    singleton = x == y
    count = length(pairs)
    last = chart.groups[end]
    if iseven(count)
        # The next odd head is normalized relative to the preceding odd head.
        hx = _chart_last_reduced(chart, pairs, x, fx, method)
        hy = _chart_last_reduced(chart, pairs, y, fy, method)
        scaled = S.(_block_line(x, hx, y, hy))
        leading = findlast(!iszero, scaled)
        ratio = leading === nothing ? one(S) : scaled[leading]
        B = scaled ./ ratio
        resize!(B, leading === nothing ? 1 : leading)
        scale = last.scale * ratio
        coefficients = last.scale .* scaled
        reduced = (last.scale * hx, last.scale * hy)
        all(isfinite, (ratio, scale, reduced..., B..., coefficients...)) || return nothing
        last.ratio = ratio
        push!(chart.groups, _QTCFGroup{S}(B, S[], zero(S),
            (zero(S), zero(S), zero(S)), one(S), scale))
        return coefficients, reduced
    end

    # A new even block closes the pending normalized odd head.
    k = length(chart.groups)
    bx, by = _chart_head(last, x), _chart_head(last, y)
    ax, ay = _chart_a(pairs[end], x), _chart_a(pairs[end], y)
    if method isa OneDiv
        ux, vx = _chart_input_pair(chart, pairs, x, fx, k)
        uy, vy = _chart_input_pair(chart, pairs, y, fy, k)
        invx = _chart_divide(one(S), ux - bx * vx)
        invy = _chart_divide(one(S), uy - by * vy)
        gx, gy = ax * (vx * invx), ay * (vy * invy)
        cx, cy = ax * (ux * invx), ay * (uy * invy)
    else
        vx = _chart_input(chart, pairs, x, fx, k)
        vy = _chart_input(chart, pairs, y, fy, k)
        gx, gy = _chart_divide(ax, vx-bx), _chart_divide(ay, vy-by)
        cx, cy = isinf(vx) ? ax : ax * _chart_divide(vx, vx-bx),
            isinf(vy) ? ay : ay * _chart_divide(vy, vy-by)
    end
    all(isfinite, (gx, gy, cx, cy)) || return nothing
    B = last.head
    gamma = S.(_block_line(x, gx, y, gy))
    center = zero(S)
    anchor = first(pairs)[1]
    if length(B) > 1 && iszero(_chart_head(last, anchor)) && !iszero(B[2])
        # Retain a complement anchored at the original selected common zero.
        center = anchor
        c0 = _chart_a(pairs[end], center)
        pair = pairs[end]
        a2 = pair[1] == pair[2] ? zero(S) : one(S)
        dx, dy = x-center, y-center
        a1 = pair[1] == pair[2] ? one(S) : (center-pair[1])+(center-pair[2])
        if singleton
            c2 = a2
            c1 = iszero(dx) ? a1+B[2]*gx : (cx-c0)/dx-c2*dx
        else
            line = S.(_block_line(dx, (cx-c0)/dx, dy, (cy-c0)/dy))
            c1, c2 = line[1], length(line) > 1 ? line[2] : zero(S)
        end
        complement = (S(c0), S(c1), S(c2))
        g0, g1 = (c1-a1)/B[2], (c2-a2)/B[2]
        gamma = S[g0-g1*center, g1]
    else
        # Other affine heads retain their polynomial coefficients without a remote anchor.
        a = _factor_poly(pairs[end]...)
        c = _padd(a, _pmul(B, gamma))
        complement = ntuple(j -> j <= length(c) ? S(c[j]) : zero(S), 3)
    end
    coefficients = gamma ./ last.scale
    reduced = (gx / last.scale, gy / last.scale)
    all(isfinite, (center, complement..., gamma..., coefficients..., reduced...)) || return nothing
    last.center = center
    last.complement = complement
    last.gamma = gamma
    return coefficients, reduced
end

# %% _new_weight
# Return the reduced value for the proposed datum using Classic or OneDiv.
function _new_weight(method::Union{Classic,OneDiv}, r::QuadraticThiele, z_new, y_new)
    d = getfield(r, :data)
    if r.chart !== nothing
        g = last(r.chart.groups)
        if method isa OneDiv
            U, V = _chart_input_pair(r.chart, r.pairs, z_new, y_new, length(r.chart.groups))
            delta = U - _chart_head(g, z_new) * V
            a = _chart_a(last(r.pairs), z_new)
            return isempty(g.gamma) ? _chart_divide(a * V, g.scale * delta) :
                _chart_divide(g.scale * a * delta,
                    _chart_complement(g, z_new) * V - _chart_gamma(g, z_new) * U)
        end
        if isempty(g.gamma)
            u = _chart_input(r.chart, r.pairs, z_new, y_new, length(r.chart.groups))
            return _chart_divide(_chart_a(last(r.pairs), z_new),
                g.scale * (u - _chart_head(g, z_new)))
        end
        return g.scale * _chart_last_reduced(r.chart, r.pairs, z_new, y_new)
    end
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
        u = iszero(D) ? convert(typeof(D), iszero(N) ? NaN : Inf) : N / D
    end
    return u
end

# %% _block_line
# Interpolate the two reduced values by a constant or affine polynomial.
function _block_line(z1, h1, z2, h2)
    (!isfinite(h1) || !isfinite(h2)) && return [oftype(h1, NaN)]
    _same(z1, z2) && return [h1]

    slope = (h2 - h1) / (z2 - z1)
    midpoint = (z1 + z2) / 2
    midpoint_value = (h1 + h2) / 2
    b = [midpoint_value - slope * midpoint, slope]

    return all(isfinite, b) ? b : eltype(b)[NaN]
end

# %% _factor_poly
# Return linear or quadratic numerator coefficients for the supplied node pair.
function _factor_poly(z1, z2)
    S = promote_type(typeof(float(z1)), typeof(float(z2)))
    singleton = _same(z1, z2)
    z1, z2 = S(z1), S(z2)
    return singleton ? S[-z1, 1] : S[z1*z2, -(z1 + z2), 1]
end

# %% [markdown]
# ## 5. Symmetry wrappers and mesh/selection helpers
#
# `Symmetric` and `symmetric` wrap point and value maps for the adapters below.
# `_clean_param`, `_reflection_param`, and `_mesh_params` locate nodes and build
# the continuum mesh.
# `_same` compares node coordinates exactly.
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
    evaluate!(rbuf, r, z, abuf, bbuf)
    @inbounds for i in 1:m
        isnan(rbuf[i]) && (rbuf[i] = evaluate(r, z[i], eval_method))
        errors[i] = abs(y[i] - rbuf[i])
    end
    return errors
end

# %% approximate — continuum sampling
"""
    approximate(f, curve, QTCF(); fixed_points, initial_node, float_type, tol, max_degree=100, kwargs...)

Use adaptive interior samples. Startup uses the parameter midpoint unless
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
    initial_refinement::Int=15, parameter_reflection::Union{Nothing,Function}=nothing)
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

    while true
        empty!(active)
        for j in 2:num_ref+1, i in 1:length(path.next)
            index = (j - 1) * rows + i
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
            [k for k in eachindex(active) if τ[active[k]] != reflection(τ[active[k]])] :
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
            reason, best = quitting_check(history, stagnation, tol, fmax, max_iter, policy)
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
The check set includes selected nodes. Every completed block is recorded;
return the lowest-error admissible record. Budgets count blocks.
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
    best = best_acceptable(history, policy)
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
    chart = _convert_chart(S, r.chart)
    if chart !== nothing
        g = first(chart.groups)
        g.head[1] += s
        gamma1 = length(g.gamma) > 1 ? g.gamma[2] : zero(S)
        g.complement = (g.complement[1] + s*_chart_gamma(g, g.center),
            g.complement[2] + s*gamma1, g.complement[3])
    end
    return _qtcf_from_blocks(Tuple{S,S}.(d.pairs), b,
        [S.(p) for p in d.apolys], S.(d.nodes), S.(d.vals .+ s), chart)
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
    chart = _convert_chart(S, r.chart)
    if chart !== nothing
        g = first(chart.groups)
        g.head .*= s
        g.gamma ./= s
        g.ratio *= s
        for g in Iterators.drop(chart.groups, 1)
            g.scale *= s
        end
    end
    return _qtcf_from_blocks(Tuple{S,S}.(d.pairs), b,
        [S.(p) for p in d.apolys], S.(d.nodes), S.(s .* d.vals), chart)
end
