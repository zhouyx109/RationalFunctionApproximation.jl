# -*- coding: utf-8 -*-
"""
    Thiele (type)

Thiele continued-fraction representation of a rational function. Aliased as `TCF`.

# Fields
- `nodes`: interpolation nodes of the rational function
- `values`: values at the nodes
- `weights`: continued-fraction (inverse-difference) coefficients
"""
struct Thiele{T,S} <: AbstractRationalInterpolant{T,S}
    nodes::Vector{S}
    values::Vector{S}
    weights::Vector{S}
    function Thiele{T}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S}
        ) where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight)
        new{T,S}(node, value, weight)
    end
    function Thiele{T,S}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S}
        ) where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight)
        new{T,S}(node, value, weight)
    end
end

# Convert the numeric type:
function Base.convert(::Type{F}, r::Thiele{T,S}) where {F<:AbstractFloat,T,S<:Real}
    return Thiele{F,F}( F.(nodes(r)), F.(values(r)), F.(weights(r)))
end

function Base.convert(::Type{F}, r::Thiele{T,S}) where {F<:AbstractFloat,T,S<:Complex}
    CF = complex(F)
    return Thiele{F,CF}( CF.(nodes(r)), CF.(values(r)), CF.(weights(r)))
end

# convenience accessors amd overloads
nodes(r::Thiele) = r.nodes
Base.values(r::Thiele) = r.values
weights(r::Thiele) = r.weights
function degrees(r::Thiele)
    n = length(r.nodes) - 1
    q, s = divrem(n, 2)
    return q + s, q
end
degree(r::Thiele) = div(length(r.nodes) - 1, 2)

Base.isreal(r::Thiele) = isreal(r.nodes) && isreal(r.weights)

function Base.copy(r::Thiele)
    return Thiele(copy(r.nodes), copy(r.values), copy(r.weights))
end

const TCF = Thiele    # alias

"""
    ThieleMethod

Abstract supertype for the recurrence strategies a [`Thiele`](@ref) can use, namely
[`OneDiv`](@ref) and [`Classic`](@ref). Pass an instance to `evaluate(r, z, method)` or
`add_node!(r, z, y, method)` for a one-off choice, or set the global default with
[`set_eval_method`](@ref) / [`set_weight_method`](@ref).

These let external or test code A/B the two algorithms. The default selectors return
singletons, so they constant-fold and the production paths stay static, allocation-free calls.
"""
abstract type ThieleMethod end

"""
    OneDiv() <: ThieleMethod

Evaluate a [`Thiele`](@ref) via the numerator/denominator pair recurrence, deferring the
single division to the end. Faster and the default for point evaluation; falls back to
[`Classic`](@ref) if the denominator underflows to zero.
"""
struct OneDiv  <: ThieleMethod end

"""
    Classic() <: ThieleMethod

Evaluate a [`Thiele`](@ref) by the textbook backward continued-fraction recurrence (one
division per node). Slower than [`OneDiv`](@ref) but the default for new-weight computation in
`add_node!`.
"""
struct Classic <: ThieleMethod end

default_eval_method()   = OneDiv()    # point evaluation
default_weight_method() = Classic()   # new-weight computation in `add_node!`

"""
    set_eval_method(m::ThieleMethod)

Globally select the recurrence used for Thiele point evaluation (`OneDiv()` or `Classic()`).
Intended for testing; each call recompiles the affected path. See also
[`set_weight_method`](@ref).
"""
set_eval_method(m::ThieleMethod)   = (@eval default_eval_method()   = $m; m)

"""
    set_weight_method(m::ThieleMethod)

Globally select the recurrence used for Thiele new-weight computation in `add_node!`
(`OneDiv()` or `Classic()`). Intended for testing; each call recompiles the affected path.
See also [`set_eval_method`](@ref).
"""
set_weight_method(m::ThieleMethod) = (@eval default_weight_method() = $m; m)

# Empty instance, used only as a method selector in `approximate`.
Thiele() = Thiele(Float64[], Float64[], Float64[])

# Evaluation at a point
function evaluate(r::Thiele, z::Number, method::ThieleMethod=default_eval_method())
    return if isinf(z)
        if isodd(length(r.nodes))
            sum(r.weights[1:2:end])
        else
            Inf
        end
    elseif isnan(z)
        NaN
    else
        _evaluate_finite_thiele(r, z, method)
    end
end

function _evaluate_finite_thiele(r, z, method)
    value = _evaluate(r, z, method)
    if !isfinite(value)
        cancelled = _evaluate_common_zero(r, z)
        cancelled === nothing || return cancelled
    end
    return value
end

# faster for array-valued evaluation
function evaluate(r::Thiele, z::AbstractArray{<:Number})
    t = similar(z, promote_type(eltype(r), eltype(z)))
    return _evaluate_thiele_samples!(t, r, z, similar(t), similar(t), default_eval_method())
end

function evaluate!(t::AbstractArray, r::Thiele, z::AbstractArray{<:Number})
    return evaluate!(t, r, z, similar(t), similar(t))
end

# In-place evaluation with caller-supplied scratch `a`, `b` (same shape as `t`),
# so a hot loop can reuse them instead of allocating on every call.
function evaluate!(t::AbstractArray, r::Thiele, z::AbstractArray{<:Number},
    a::AbstractArray, b::AbstractArray; fallback=false)
    # use 3-term pair recurrence to avoid division until the end
    retry = fallback ? fill!(similar(t, Bool), false) : nothing
    n = length(r.weights)
    if n == 1
        t .= r.weights[1]
    else
        a .= r.weights[n]
        b .= z .- r.nodes[n-1]
        _mark_thiele_underflow!(retry, a, b)
        @inbounds for k in n-1:-1:2
            b .= r.weights[k] .* a .+ b
            a .*= z .- r.nodes[k-1]
            a, b = b, a
            _mark_thiele_underflow!(retry, a, b)
        end
        t .= r.weights[1] .+ b ./ a
    end
    return _finish_thiele_array!(t, r, z, retry)
end

_mark_thiele_underflow!(::Nothing, a, b) = nothing

function _mark_thiele_underflow!(retry, a, b)
    @inbounds for i in eachindex(retry, a, b)
        retry[i] |= guard(a[i]) || guard(b[i])
    end
    return nothing
end

function _finish_thiele_array!(t, r, z, retry)
    # Retry selected common zeros and exceptional inputs through scalar evaluation.
    @inbounds for i in eachindex(t, z)
        if isinf(z[i]) || (retry !== nothing && (retry[i] || isnan(z[i]) || !isfinite(t[i]))) ||
            (!isfinite(t[i]) && any(==(z[i]), nodes(r)))
            t[i] = evaluate(r, z[i], OneDiv())
        end
    end
    return t
end

_evaluate_thiele_samples!(t, r, z, a, b, ::OneDiv) = evaluate!(t, r, z, a, b; fallback=true)

function _evaluate_thiele_samples!(t, r, z, a, b, method::Classic)
    @inbounds for i in eachindex(t, z)
        t[i] = evaluate(r, z[i], method)
    end
    return t
end

function _evaluate_classic(r::Thiele, z::Number)
    n = length(r.nodes)
    u = last(r.weights)
    @inbounds for k in n-1:-1:1
        u = r.weights[k] + (z - r.nodes[k]) / u
    end
    return u
end

function _evaluate_onediv(r::AbstractRationalInterpolant, z::Number)
    numer, denom = _evaluate_numden(r, z)
    return if iszero(denom)
        @debug "Evaluation produced a division by zero at " z
        _evaluate_classic(r, z)  # fallback to slow evaluation
    else
        numer / denom
    end
end

_evaluate(r::AbstractRationalInterpolant, z, ::OneDiv)  = _evaluate_onediv(r, z)
_evaluate(r::AbstractRationalInterpolant, z, ::Classic) = _evaluate_classic(r, z)

guard(z) = !iszero(z) && abs(real(z)) < 1e-100 && abs(imag(z)) < 1e-100
function _evaluate_numden(r::Thiele, z::Number; underflow_guard=true)
    # use 3-term pair recurrence to avoid division until the end
    @assert isfinite(z)
    n = length(r.weights)
    return if n == 1
        r.weights[1], 1
    else
        a = r.weights[n]
        b = z - r.nodes[n-1]
        @inbounds for k in n-1:-1:2
            t = r.weights[k] * a + b
            b = a * (z - r.nodes[k-1])
            a = t
            if underflow_guard && guard(a)    # prevent underflow
                @debug "scaling to avoid underflow at " z
                a *= 1e100
                b *= 1e100
            end
        end
        r.weights[1] * a + b, a
    end
end

# A scalar node is a block with b_k=w_k and a_k(t)=t-z_k.
function _evaluate_common_zero(r::Thiele, z)
    j = findfirst(==(z), r.nodes)
    m = length(r.weights)
    (j === nothing || j == m) && return nothing
    S = promote_type(eltype(r), typeof(z))
    bjet = k -> (r.weights[k], zero(S))
    ajet = k -> (z - r.nodes[k], one(S))
    return _cancel_thiele_factor(S, j, m, bjet, ajet)
end

# Evaluate after dividing the two lower-tail polynomials by t-z at block j.
# Each jet supplies the value and first derivative of its block at z.
function _cancel_thiele_factor(::Type{S}, j, m, bjet, ajet) where {S}
    p0, p1 = bjet(m)
    p0, p1 = S(p0), S(p1)
    q0, q1 = one(S), zero(S)
    for k in m-1:-1:j+1
        b0, b1 = bjet(k)
        a0, a1 = ajet(k)
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

    bj, _ = bjet(j)
    _, factor = ajet(j)
    P, Q = bj * p1 + factor * q0, p1
    for k in j-1:-1:1
        bk, _ = bjet(k)
        ak, _ = ajet(k)
        P, Q = bk * P + ak * Q, P
        if k > 1 && guard(P)
            P *= 1e100
            Q *= 1e100
        end
    end
    return all(isfinite, (P, Q)) && !(iszero(P) && iszero(Q)) ? P / Q : nothing
end

function _thiele_jets!(ζ, r::Thiele, order, A, B)
    w, z = weights(r), nodes(r)
    n = length(w)
    if n <= 1
        A[1] = n == 0 ? zero(eltype(A)) : w[1]
        B[1] = one(eltype(B))
        A[2:order+1] .= 0; B[2:order+1] .= 0
        return nothing
    end
    A[1] = w[n]
    A[2:order+1] .= 0
    B[1] = ζ - z[n-1]
    order > 0 && (B[2] = 1)
    B[3:order+1] .= 0
    @inbounds for k in n-1:-1:2
        d = ζ - z[k-1]
        for m in order:-1:0
            t = w[k] * A[m+1] + B[m+1]
            B[m+1] = d * A[m+1] + (m > 0 ? m * A[m] : 0)
            A[m+1] = t
        end
        _guard_thiele_jets!(view(A, 1:order+1), view(B, 1:order+1))
    end
    for j in 1:order+1
        A[j], B[j] = w[1] * A[j] + B[j], A[j]
    end
    return nothing
end

function _guard_thiele_jets!(A, B)
    if any(guard, A)
        A .*= 1e50
        B .*= 1e50
    end
    return nothing
end

function _thiele_quotient_derivatives!(z, r, order, A, B, vals)
    _thiele_jets!(z, r, order, A, B)
    # A common simple factor shifts both Taylor series by one term.
    if iszero(A[1]) && iszero(B[1]) && any(==(z), nodes(r))
        _thiele_jets!(z, r, order+1, A, B)
        if all(isfinite, (A[2], B[2])) && !(iszero(A[2]) && iszero(B[2]))
            for q in 0:order
                A[q+1], B[q+1] = A[q+2]/(q+1), B[q+2]/(q+1)
            end
        end
    end
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

function _thiele_derivative(r, orders)
    isempty(orders) && throw(ArgumentError("at least one derivative order is required"))
    all(>=(0), orders) || throw(ArgumentError("derivative orders must be nonnegative"))
    order, index = maximum(orders), 1 .+ orders
    A = Vector{complex(eltype(r))}(undef, order + 2)
    B = similar(A)
    vals = similar(A, order + 1)
    return function(z)
        _thiele_quotient_derivatives!(z, r, order, A, B, vals)
        result = vals[index]
        return isreal(z) && isreal(r) ? real(result) : result
    end
end

function derivative(r::Thiele, order::Integer=1)
    f = derivative(r, [order])
    return z -> only(f(z))
end

derivative(r::Thiele, orders::AbstractVector{<:Integer}) = _thiele_derivative(r, orders)


# Return commonly normalized numerator/denominator jets for residue evaluation.
function _normalize_residue_jets(z, a, b, aʹ, bʹ)
    scale = max(one(abs(z)), abs(a), abs(b), abs(aʹ), abs(bʹ))
    return a/scale, b/scale, aʹ/scale, bʹ/scale
end

function _evaluate_numden_derivs(r::Thiele, z::Number)
    n = length(r.weights)
    return if n == 1
        r.weights[1], 1, 0, 0
    else
        a = r.weights[n]
        b = z - r.nodes[n-1]
        aʹ = 0
        bʹ = 1
        a, b, aʹ, bʹ = _normalize_residue_jets(z, a, b, aʹ, bʹ)
        @inbounds for k in n-1:-1:2
            d = z - r.nodes[k-1]
            t = r.weights[k] * aʹ + bʹ
            bʹ = a + d * aʹ
            aʹ = t
            t = r.weights[k] * a + b
            b = a * d
            a = t
            a, b, aʹ, bʹ = _normalize_residue_jets(z, a, b, aʹ, bʹ)
        end
        r.weights[1] * a + b, a, r.weights[1] * aʹ + bʹ, aʹ
    end
end

# Apply the same power-of-two row and column factors to both pencil matrices.
function _equilibrate_thiele_pencil!(D, C)
    T = real_type(promote_type(eltype(D), eltype(C)))
    for _ in 1:8
        changed = false
        for dim in (1, 2), i in axes(D, dim)
            d, c = selectdim(D, dim, i), selectdim(C, dim, i)
            magnitude = max(maximum(abs, d), maximum(abs, c))
            isfinite(magnitude) && !iszero(magnitude) || continue
            exponent = -round(Int, log2(magnitude) / 2)
            iszero(exponent) && continue
            factor = ldexp(one(T), exponent)
            if exponent < 0
                # Leave a row or column intact if scaling would round a subnormal.
                inverse = ldexp(one(T), -exponent)
                all(x -> x * factor * inverse == x, d) || continue
                all(x -> x * factor * inverse == x, c) || continue
            end
            d .*= factor
            c .*= factor
            changed = true
        end
        changed || break
    end
    return nothing
end

# The pencil builders own these matrices, so they can be scaled in place.
function _scaled_thiele_eigvals!(D, C; finite_only=true)
    _equilibrate_thiele_pencil!(D, C)
    return _thiele_eigvals(D, C; finite_only)
end

function _thiele_eigvals(D, C; finite_only=true)
    spectrum = try
        eigvals(D, C)
    catch
        _, _, _, _, ⍺, β = schur(complex(D), complex(C))
        ⍺ ./ β
    end
    return finite_only ? filter(isfinite, spectrum) : spectrum
end

function _pole_newton_step(r, z)
    _, q, _, derivative = _evaluate_numden_derivs(r, z)
    return q / derivative
end

# Refine each finite eigenvalue against the fitted denominator in its own basis.
function _refine_thiele_pencil_values!(r, spectrum)
    for i in eachindex(spectrum)
        z = spectrum[i]
        isfinite(z) || continue
        step = _pole_newton_step(r, z)
        for _ in 1:8
            isfinite(step) || break
            candidate = z - step
            isfinite(candidate) && candidate != z || break
            next_step = _pole_newton_step(r, candidate)
            isfinite(next_step) && abs(next_step) < abs(step) || break
            z, step = candidate, next_step
        end
        spectrum[i] = z
    end
    return spectrum
end

function _thiele_pencil_values(r, A, B; scale, refine, finite_only)
    spectrum = scale ? _scaled_thiele_eigvals!(A, B; finite_only) :
        _thiele_eigvals(A, B; finite_only)
    return refine ? _refine_thiele_pencil_values!(r, spectrum) : spectrum
end

# Eliminate alternating diagonal pivots and split each resulting quadratic
# off-diagonal product into two affine entries with the same product.
function _compact_thiele_pencil(r::Thiele{T,S}, eliminated) where {T,S}
    n = length(r.nodes)
    all(j -> isfinite(r.weights[j]) && !iszero(r.weights[j]), eliminated) || return nothing
    first_kept = 5 - first(eliminated)
    kept = first_kept:2:n
    m = length(kept)
    D, C = zeros(S, m, m), zeros(S, m, m)
    for (i, j) in enumerate(kept)
        D[i, i] = r.weights[j]
        if j > 2
            pivot = r.weights[j-1]
            D[i, i] -= r.nodes[j-1] / pivot
            C[i, i] -= inv(pivot)
        end
        if j < n
            pivot = r.weights[j+1]
            D[i, i] -= r.nodes[j] / pivot
            C[i, i] -= inv(pivot)
            if i < m
                D[i, i+1] = -r.nodes[j] / pivot
                D[i+1, i] = -r.nodes[j+1] / pivot
                C[i, i+1] = C[i+1, i] = -inv(pivot)
            end
        end
    end
    return all(isfinite, D) && all(isfinite, C) ? (D, C) : nothing
end

"""
    poles(r::Thiele; compact=false, finite_only=true, scale=false, refine=false)

Return the finite roots of the unreduced denominator from the full scalar pencil,
solved as a generalized eigenvalue problem with no scaling and no refinement.
The keywords are opt-in departures from that default. With `compact=true`,
try eliminating alternating rows and columns before solving; if neither set of
pivots yields a finite compact pencil, use the full pencil.
Set `scale=true` for joint power-of-two row and column scaling before the solve, or
`refine=true` for up to eight guarded Newton steps against the fitted denominator
at each finite eigenvalue.
Set `finite_only=false` to return every eigensolver value, including Inf and NaN.
"""
function poles(r::Thiele{T,S}; compact=false, finite_only=true,
    scale::Bool=false, refine::Bool=false) where {T,S}
    n = length(r.nodes)
    (n < 2 || ((finite_only || compact) && n == 2)) && return S[]
    if compact
        pencil = _compact_thiele_pencil(r, 2:2:n)
        if isnothing(pencil)
            pencil = _compact_thiele_pencil(r, 3:2:n)
        end
        if !isnothing(pencil)
            return _thiele_pencil_values(r, pencil...; scale, refine, finite_only)
        end
    end
    C = diagm(1 => -ones(T, n-2))
    xi = view(r.nodes, 2:n-1)
    D = diagm(0 => view(r.weights, 2:n), -1 => -ones(T, n-2), 1 => -xi)
    return _thiele_pencil_values(r, D, C; scale, refine, finite_only)
end

# Finite numerator roots from the full scalar pencil; `scale=true` opts in to pencil scaling.
function roots(r::Thiele{S,T}; scale::Bool=false) where {S,T}
    n = length(r.nodes)
    n == 0 && return T[]
    C = diagm(1 => -ones(T, n-1))
    xi = view(r.nodes, 1:n-1)
    D = diagm(0 => r.weights, -1 => -ones(T, n-1), 1 => -xi)
    return scale ? _scaled_thiele_eigvals!(D, C) : _thiele_eigvals(D, C)
end

# Compute TCF poles and pass them to the shared residue loop.
residues(r::Thiele; check_finite=false, radius=nothing, compact=false,
    scale::Bool=false, refine::Bool=false) =
    residues(r, poles(r; compact, scale, refine); check_finite, radius)

# Extract residues at supplied poles; a callback can reuse numerator/denominator workspace.
function residues(r::AbstractRationalInterpolant, ζ;
    check_finite=false, radius=nothing, numden=_evaluate_numden_derivs)
    res = similar(ζ, complex(eltype(ζ)))
    for i in eachindex(ζ)
        # is it a simple pole?
        p, q, pʹ, qʹ = numden(r, ζ[i])
        candidate = iszero(qʹ) ? oftype(p, NaN) : p / qʹ
        if !iszero(qʹ) && (!check_finite || isfinite(candidate))
            # simple pole
            res[i] = candidate
        else
            @debug "Fallback for residue at pole " ζ[i]
            # not a simple pole; use fallback contour integral
            res[i] = radius === nothing ? Res(r, ζ[i]; avoid=ζ) :
                Res(r, ζ[i]; avoid=ζ, radius=radius(ζ[i]))
        end
    end
    return ζ, res
end

function Thiele(nodes::AbstractVector, values::AbstractVector, weights::AbstractVector)
    if isempty(nodes) && eltype(nodes) == Any
        nodes = values = weights = Float64[]
    end
    nodes, values, weights = promote(float(nodes), float(values), float(weights))
    T = eltype(values)
    return Thiele{real_type(T)}(nodes, values, weights)
end

function Thiele(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(DimensionMismatch("node and value counts must agree"))
    idx = axes(x, 1)
    d = collect(float(promote_type(eltype(x), eltype(y))), y)
    for (len, i) in enumerate(idx)
        for k in first(idx, len-1)
            d[i] = (x[i] - x[k]) / (d[i] - d[k])
        end
    end
    return Thiele(collect(x), y, d)
end

# update in-place the nodes and weights
function add_node!(r::Thiele, z_new, y_new, method::ThieleMethod=default_weight_method())
    isfinite(z_new) && isfinite(y_new) || throw(ArgumentError("Interpolation data must be finite"))
    w = _new_weight(method, r.nodes, r.weights, z_new, y_new)
    isnan(w) && throw(NaNException("Adding node at $z_new caused NaN weight"))
    isfinite(w) || throw(ArgumentError("Adding node at $z_new caused an infinite weight"))
    push!(r.weights, w)
    push!(r.nodes, z_new)
    push!(r.values, y_new)
    return r
end

function _new_weight_onediv(z, w, z_new, y_new)
    a = 1
    b = y_new
    @inbounds for k in eachindex(z)
        t = -w[k] * a + b
        b = a * (z_new - z[k])
        a = t
    end
    return b / a
end

function _new_weight_classic(z, w, z_new, y_new)
    u = y_new
    @inbounds for k in eachindex(z)
        u = _classic_reduced_value(z_new - z[k], u - w[k])
    end
    return u
end

_classic_reduced_value(N, D) = iszero(D) ? convert(typeof(D), iszero(N) ? NaN : Inf) : N / D

_new_weight(::OneDiv,  z, w, z_new, y_new) = _new_weight_onediv(z, w, z_new, y_new)
_new_weight(::Classic, z, w, z_new, y_new) = _new_weight_classic(z, w, z_new, y_new)

# Per-iteration kernels for the continuum `approximate` loop below, gathered under one name.
# They are kept as separate functions so each specializes on the concrete element types of its
# array arguments: even though `DiscretizedPath.points` is parameterized, ComplexRegions' curve
# evaluation (`point(d, ·)` / `path(t)`) infers as `Any`, so `σ`, `fσ`, and everything derived
# from them (`fτ`, `rbuf`, `zbuf`) have non-inferable element types inside `approximate`. Inline
# per-element loops over them would dispatch dynamically on every element; routing through these
# barriers replaces that with one dispatch per call plus a type-stable loop.

# Evaluate r at the active test points (reusing the scratch buffers) and return the maximum
# error and the index (into `active`) of the worst point; mirrors `findmax` (NaN wins).
_eligible_thiele_pair(z, partner, selected) = isfinite(z) && isfinite(partner) &&
    z != partner && (partner + zero(partner)) ∉ selected

function _sweep!(rbuf, zbuf, abuf, bbuf, r::Thiele, τ, fτ, active,
    eval_method::ThieleMethod=default_eval_method())
    m = length(active)
    resize!(rbuf, m); resize!(zbuf, m); resize!(abuf, m); resize!(bbuf, m)
    @inbounds for i in 1:m
        zbuf[i] = τ[active[i]]
    end
    _evaluate_thiele_samples!(rbuf, r, zbuf, abuf, bbuf, eval_method)
    err_max = abs(fτ[active[1]] - rbuf[1])
    kmax = 1
    @inbounds for k in 2:m
        e = abs(fτ[active[k]] - rbuf[k])
        if isnan(e) || e > err_max
            err_max, kmax = e, k
            isnan(e) && break
        end
    end
    return err_max, kmax
end

# Sample f at the points indexed by `idx`, storing into `fτ`; return max |f| over them.
function _sweep!(fτ, f::Function, τ, idx)
    fmax = abs(zero(eltype(fτ)))
    @inbounds for i in idx
        v = f(τ[i])
        fτ[i] = v
        a = abs(v)
        (a > fmax) && (fmax = a)
    end
    return fmax
end

"""
    approximate(f, curve, TCF(); reflection=nothing, kwargs...)

Fit a scalar Thiele fraction on an adaptive interior mesh. With a reflection,
start at the reflection-fixed parameter midpoint. Insert each selected point
and its reflected partner using two scalar weights. Record only completed
pairs and stop if either weight is nonfinite. In this mode `max_iter` counts
the midpoint and completed pairs, and refinement decreases once per pair.

`value_reflection(z, fz, w)` may supply the value at the reflected point
`w`. `parameter_reflection(s)` may supply its mesh parameter. The stopping
test uses the sampled target scale unless `comparison_scale` is supplied.
`weight_method` selects scalar weight construction and `eval_method` selects
the sweep evaluator. OneDiv retries suspect array results with scalar evaluation.
For a nonfinite paired weight, `result.options.failure` records the block,
leg, node and weight. It is `nothing` otherwise.
"""
function approximate(
    f::Function, d::Union{ComplexPath,ComplexCurve}, ::Thiele;
    float_type::Type = promote_type(real_type(d), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed = true,
    max_degree = 100,
    reflection::Union{Nothing,Function}=nothing,
    max_iter::Int = reflection === nothing ? 2max_degree : max_degree,
    refinement::Int = 3,
    stagnation::Int = 5,
    initial_refinement::Int=15,
    value_reflection::Union{Nothing,Function}=nothing,
    parameter_reflection::Union{Nothing,Function}=nothing,
    comparison_scale=nothing,
    weight_method::ThieleMethod=default_weight_method(),
    eval_method::ThieleMethod=default_eval_method()
    )

    if reflection !== nothing
        return _approximate_paired_thiele(f, d, reflection;
            value_reflection, parameter_reflection, comparison_scale, tol, allowed,
            max_iter, refinement, initial_refinement, stagnation, weight_method, eval_method)
    end
    if value_reflection !== nothing || parameter_reflection !== nothing || comparison_scale !== nothing
        throw(ArgumentError("value_reflection, parameter_reflection and comparison_scale require reflection"))
    end
    if allowed == :strict
        # only allow poles off the curve
        allowed = z -> dist(z, d) > tol
    end

    num_ref = initial_refinement
    path = DiscretizedPath(d, [0, 1]; refinement=num_ref, maxpoints=max_iter * refinement)
    σ = [point(d, 0)]
    fσ = f.(σ)         # f at nodes

    # Test points have one matrix row per node (except the last). `active` holds *linear*
    # indices of the currently live test points; it need not be rectangular. We grow it
    # incrementally (each accepted node appends one row's worth of test points) and only
    # rebuild it wholesale when the path is re-discretized in the initial phase. Linear
    # integer indexing is much cheaper here than CartesianIndex indexing, which dominated
    # the per-point loops below via `to_indices`.
    # `path.points` is concretely typed at runtime, but inference here can't prove it (its
    # element type flows from `point(d, ·)`, which infers as `Any`), so the per-point work
    # goes through the kernel barriers above rather than relying on the type in this scope.
    τ = path.points
    R = size(τ, 1)        # rows; linear index of (i, j) is (j-1)*R + i
    active = Int[]
    sizehint!(active, max_iter * refinement)
    for j in 2:num_ref+1
        push!(active, (j - 1) * R + 1)       # row 1, columns 2:num_ref+1
    end
    fτ = Matrix{eltype(fσ)}(undef, size(τ))        # f at test points
    fmax = _sweep!(fτ, f, τ, active)        # sample f; fmax is the scale of f

    # Reusable scratch, aligned with `active`: the gathered test points `zbuf`, r there
    # (`rbuf`), and the two work arrays the recurrence in `evaluate!` would otherwise
    # allocate on every call. Gathering into a dense `zbuf` lets `evaluate!` run on
    # contiguous vectors. Use the natural recurrence type (values promoted with the
    # nodes/points) so real problems stay real-valued instead of paying for complex
    # arithmetic, while complex domains (e.g. an imaginary interval) still get a wide
    # enough buffer.
    S = promote_type(eltype(fσ), eltype(σ))
    zbuf = Vector{eltype(σ)}()
    rbuf = S[]
    abuf = S[]
    bbuf = S[]

    # Initialize rational approximation
    r = Thiele(σ, fσ)
    history = [IterationRecord(r, NaN, missing)]

    # Main iteration
    n = 1       # iteration counter
    stop = nothing
    while true
        err_max, idx_max = _sweep!(rbuf, zbuf, abuf, bbuf, r, τ, fτ, active, eval_method)
        history[n].error = err_max

        reason, best = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if reason !== :iterating
            if reason !== :converged
                @info("Stopping with estimated error $(round(history[best].error, sigdigits=4)) after $n iterations")
                r = history[best].interpolant
            end
            stop = ConvergenceStatus(reason, best, history)
            break
        end

        # Add node to approximant
        idx_new = active[idx_max]      # location of worst test point
        try
            add_node!(r, τ[idx_new], fτ[idx_new], weight_method)
            push!(history, IterationRecord(r, NaN, missing))
        catch e
            e isa Union{ArgumentError,NaNException} || rethrow()
            best = best_acceptable(history, allowed)
            r = history[best].interpolant
            reason = e isa NaNException ? :nan_weight : :node_failure
            stop = ConvergenceStatus(reason, best, history)
            @info("Adding node failed; stopping with estimated error $(round(history[best].error, sigdigits=4))")
            @debug("Error $e")
            break
        end

        # Add node to path (needs the (row, column) form of the worst point)
        local idx_new_test
        try
            idx_new_test = add_node!(path, CartesianIndices(τ)[idx_new])
        catch
            # look for the best acceptable case
            best = best_acceptable(history, allowed)
            r = history[best].interpolant
            stop = ConvergenceStatus(:refinement, best, history)
            @info("Maximum path refinement exceeded; stopping with estimated error $(round(history[best].error, sigdigits=4))")
            break
        end

        n += 1
        # In the initial phase, we throw out the old test points.
        if num_ref > refinement    # initial phase
            num_ref = max(refinement, num_ref - 1)    # gradually decrease refinement level
            s = first(collect(path))
            reset!(path, s; refinement=num_ref)
            empty!(active)
            for j in 2:num_ref+1, i in 1:n
                push!(active, (j - 1) * R + i)
            end
            fmax = _sweep!(fτ, f, τ, active)
        else
            # Evaluate f at the new test points only, and extend `active` by the new row.
            _sweep!(fτ, f, τ, idx_new_test)
            new_row = idx_new_test[end][1]
            for j in 2:num_ref+1
                push!(active, (j - 1) * R + new_row)
            end
        end
    end
    return ContinuumApproximation(f, d, r, allowed, path, history, stop)
end

function _approximate_paired_thiele(f, d, reflection;
    value_reflection, parameter_reflection, comparison_scale, tol, allowed,
    max_iter, refinement, initial_refinement, stagnation, weight_method, eval_method)
    refinement >= 1 || throw(ArgumentError("refinement must be positive"))
    max_iter >= 1 || throw(ArgumentError("max_iter must be positive"))
    T = real_type(d)
    terminal = T(length(d))
    midpoint = terminal/2
    z0 = point(d, midpoint)
    reflection(z0) == z0 || throw(ArgumentError("The parameter midpoint must be fixed by reflection"))
    y0 = f(z0)
    isfinite(y0) || throw(ArgumentError("The midpoint value must be finite"))
    r = Thiele([z0], [y0], [y0])
    S = eltype(r)
    policy = _qtcf_allowed(allowed, d, tol)
    sampler = _QTCFCurveSampler(d, f, parameter_reflection)
    path = DiscretizedPath(d, T[0, terminal]; refinement=max(refinement, initial_refinement),
        maxpoints=max(2max_iter+2, max_iter*refinement))
    selected_params = T[midpoint]
    selected_values = S[y0]
    reset!(path, _mesh_params(d, selected_params, T); refinement=path.nref)
    τ, num_ref = path.points, path.nref
    rows = size(τ, 1)
    active = Int[]
    eligible = Int[]
    fτ = Matrix{S}(undef, size(τ))
    fmax = zero(T)
    ztest, ytest, rbuf, abuf, bbuf = S[], S[], S[], S[], S[]
    errors = real_type(S)[]
    history = [IterationRecord(r, NaN, missing)]
    reason = :max_degree
    failure = nothing

    while true
        empty!(active)
        for j in 2:num_ref+1, i in eachindex(path.next)
            index = (j - 1) * rows + i
            0 < path.params[index] < terminal && τ[index] ∉ nodes(r) && push!(active, index)
        end
        _sweep!(fτ, f, τ, active)
        empty!(ztest); empty!(ytest)
        append!(ztest, view(τ, active))
        append!(ytest, view(fτ, active))
        for i in eachindex(selected_params)
            push!(ztest, nodes(r)[i])
            push!(ytest, selected_values[i])
        end
        for buffer in (rbuf, abuf, bbuf)
            resize!(buffer, length(ztest))
        end
        resize!(errors, length(ztest))
        _evaluate_thiele_samples!(rbuf, r, ztest, abuf, bbuf, eval_method)
        for i in eachindex(errors)
            errors[i] = abs(ytest[i] - rbuf[i])
        end
        history[end].error = maximum(errors)
        fmax = max(fmax, maximum(abs, ytest))
        reason, _ = quitting_check(history, stagnation, tol, something(comparison_scale, fmax),
            max_iter, policy)
        reason === :iterating || break

        empty!(eligible)
        for k in eachindex(active)
            z = τ[active[k]]
            partner = reflection(z)
            _eligible_thiele_pair(z, partner, nodes(r)) && push!(eligible, k)
        end
        if isempty(eligible)
            reason = :exhausted
            break
        end
        _, worst = findmax(view(errors, eligible))
        index = active[eligible[worst]]
        s, z, y = path.params[index], τ[index], fτ[index]
        block = _qtcf_block(sampler, s, z, y, reflection, value_reflection)
        if block === nothing || !all(isfinite, block.values)
            reason = :nonfinite_target
            break
        end
        trial = copy(r)
        for (leg, (node, value)) in enumerate(zip(block.nodes, block.values))
            weight = _new_weight(weight_method, trial.nodes, trial.weights, node, value)
            if !isfinite(weight)
                failure = (; block=length(history)+1, leg, node, weight)
                break
            end
            push!(trial.nodes, node)
            push!(trial.values, value)
            push!(trial.weights, weight)
        end
        if failure !== nothing
            reason = :nonfinite_weight
            break
        end
        r = trial
        push!(history, IterationRecord(r, NaN, missing))
        append!(selected_params, block.ids)
        append!(selected_values, block.raw)
        try
            num_ref = max(refinement, num_ref-1)
            reset!(path, _mesh_params(d, selected_params, T); refinement=num_ref)
        catch error
            error isa Union{ArgumentError,KeyError,MaxRefinementException} || rethrow()
            reason = :refinement
            break
        end
    end
    best = best_acceptable(history, policy)
    r = copy(history[best].interpolant)
    stop = ConvergenceStatus(reason, best, history)
    return ContinuumApproximation(f, d, r, policy, path, history, stop; options=(;failure))
end

function approximate(
    y::AbstractVector{T}, z::AbstractVector{S}, ::Thiele;
    float_type::Type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::AbstractFloat = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = min(length(y), 200),
    stagnation::Int = 5,
    ) where {T<:Number,S<:Number}

    length(y) == length(z) || throw(DimensionMismatch("node and value counts must agree"))
    isempty(z) && throw(ArgumentError("At least one sample is required"))
    fmax = maximum(abs, y)     # scale of f
    _, idx_min = findmin(abs, y)
    r = Thiele([z[idx_min]], [y[idx_min]])
    y_test = copy(collect(y))
    z_test = copy(collect(z))
    deleteat!(y_test, idx_min)
    deleteat!(z_test, idx_min)
    r_test = similar(y_test, eltype(r))
    abuf, bbuf = similar(r_test), similar(r_test)

    history = [IterationRecord(r, NaN, missing)]
    n = 1    # iteration counter
    stop = nothing
    while true
        if isempty(z_test)
            history[n].error = maximum(abs(r(z[i]) - y[i]) for i in eachindex(y, z))
            reason, best = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
            if reason !== :converged
                reason, best = :exhausted, best_acceptable(history, allowed)
                r = history[best].interpolant
            end
            stop = ConvergenceStatus(reason, best, history)
            break
        end
        resize!(abuf, length(r_test)); resize!(bbuf, length(r_test))
        _evaluate_thiele_samples!(r_test, r, z_test, abuf, bbuf, default_eval_method())
        r_test .-= y_test
        err_max, idx_max = findmax(abs(e) for e in r_test)
        history[n].error = err_max

        reason, best = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if reason !== :iterating
            if reason !== :converged
                @info("Stopping with estimated error $(round(history[best].error, sigdigits=4)) after $n iterations")
                r = history[best].interpolant
            end
            stop = ConvergenceStatus(reason, best, history)
            break
        end

        # Add new node:
        try
            add_node!(r, z_test[idx_max], y_test[idx_max])
            push!(history, IterationRecord(r, NaN, missing))
            deleteat!(y_test, idx_max)
            deleteat!(z_test, idx_max)
            deleteat!(r_test, idx_max)
        catch(e)
            # look for the best acceptable case
            best = best_acceptable(history, allowed)
            r = history[best].interpolant
            stop = ConvergenceStatus(:node_failure, best, history)
            @info("Adding node failed; stopping with estimated error $(round(history[best].error, sigdigits=4))")
            @debug("Error $e")
            break
        end
        n += 1
    end
    return DiscreteApproximation(y, z, r, _remaining_mask(r, z), allowed, history, stop)
end

# Operations with scalars that can be done quickly.

function Base.:+(r::Thiele, s::Number)
    w = copy(r.weights)
    w = [w[1] + s; w[2:end]]
    return Thiele(r.nodes, r.values .+ s, w)
end

function Base.:-(r::Thiele)
    return Thiele(r.nodes, -r.values, -r.weights)
end

function Base.:*(r::Thiele, s::Number)
    return if iszero(s)
        zer = zero(eltype(r.values))
        Thiele(r.nodes[1:1], [zer], [zer])
    else
        y = s * r.values
        w = Vector{eltype(y)}(undef, length(r.weights))
        w[1:2:end] .= r.weights[1:2:end] * s
        w[2:2:end] .= r.weights[2:2:end] / s
        Thiele(r.nodes, y, w)
    end
end
