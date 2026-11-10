# Construct QTCF blocks in normalized coordinates for this comparison.
import RationalFunctionApproximation as RFA
using LinearAlgebra: eigvals, schur
using RationalFunctionApproximation: Classic, QuadraticThiele, nodes

struct NormalizedQTCF <: RFA.ThieleMethod
    trace::Union{Nothing,Vector{NamedTuple}}
end
NormalizedQTCF() = NormalizedQTCF(nothing)

RFA._new_weight(::NormalizedQTCF, z, w, zn, yn) =
    RFA._new_weight(Classic(), z, w, zn, yn)

# Replay an inverse block on a ratio represented by two coordinates.
function normalized_coordinates(r, z, fz)
    n, d = fz, one(fz)
    for (a, b, pair) in zip(r.apolys, r.bpolys, r.pairs)
        A = RFA._evaluate_block(a, pair, z; numerator=true)
        B = RFA._evaluate_block(b, pair, z)
        n, d = A * d, n - B * d
        magnitude = max(abs(n), abs(d))
        if isfinite(magnitude) && !iszero(magnitude) &&
            !(0x1p-100 <= magnitude <= 0x1p100)
            scale = exp2(clamp(exponent(magnitude), -1022, 1023))
            n, d = n / scale, d / scale
        end
    end
    return n, d
end

function normalized_values(r, x, w, fx, fw)
    ratios = (normalized_coordinates(r, x, fx), normalized_coordinates(r, w, fw))
    all(pair -> all(isfinite, pair) && !iszero(last(pair)), ratios) ||
        throw(ArgumentError("zero or nonfinite projective denominator"))
    # Choose a common factor before taking either ratio.
    log_magnitude = maximum(log2(abs(n)) - log2(abs(d)) for (n, d) in ratios)
    factor = isinf(log_magnitude) ? 1.0 : exp2(clamp(-ceil(log_magnitude), -1000, 1000))
    hx, hw = map(pair -> (factor * first(pair)) / last(pair), ratios)
    return RFA._block_line(x, hx, w, hw), factor, (hx, hw)
end

function RFA.add_node!(r::QuadraticThiele, z::AbstractVector, y::AbstractVector,
    method::NormalizedQTCF)
    length(z) == length(y) == 2 ||
        throw(ArgumentError("normalized QTCF requires a pair after midpoint startup"))
    S = eltype(r)
    x, w = S.(z)
    fx, fw = S.(y)
    all(isfinite, (x, w, fx, fw)) || throw(ArgumentError("nonfinite pair data"))
    (x == w || x in nodes(r) || w in nodes(r)) &&
        throw(ArgumentError("normalized QTCF requires unused distinct nodes"))
    block = length(r.pairs) + 1
    try
        b, factor, ordinates = normalized_values(r, x, w, fx, fw)
        scaled = factor .* last(r.apolys)
        all(isfinite, (b..., scaled..., ordinates...)) && !iszero(last(scaled)) ||
            throw(ArgumentError("normalized block coefficients out of range"))
        RFA._qtcf_append_block!(r, (x, w), (fx, fw), b)
        # This numerator was unused in the previous convergent.
        r.apolys[end-1] = scaled
        if method.trace !== nothing
            push!(method.trace, (;block, reduced_max=maximum(abs, ordinates),
                intermediate_max=NaN, slope=length(b) == 2 ? abs(b[2]) : 0.0,
                factor=abs(factor), outcome="normalized", detail=""))
        end
        return r
    catch e
        if method.trace !== nothing && e isa Union{ArgumentError,RFA.NaNException}
            push!(method.trace, (;block, reduced_max=NaN, intermediate_max=NaN,
                slope=NaN, factor=NaN, outcome="failed", detail=sprint(showerror, e)))
        end
        rethrow()
    end
end

# Split each nonmonic numerator across the two off-diagonals, keeping their
# product equal to its coefficient. Residues use the same stored fraction.
function normalized_residues(r::QuadraticThiele)
    S = eltype(r)
    m = length(r.pairs) - 1
    m > 0 || return S[], S[]
    C, D = zeros(S, m, m), zeros(S, m, m)
    for j in 1:m
        b = r.bpolys[j+1]
        C[j, j] = b[1]
        length(b) == 2 && (D[j, j] = b[2])
        j == 1 && continue
        a = r.apolys[j]
        coefficient = last(a)
        left = sqrt(abs(coefficient))
        right = coefficient / left
        x, w = r.pairs[j]
        C[j-1, j], D[j-1, j] = -left*x, left
        if length(a) == 2
            C[j, j-1] = -right
        else
            C[j, j-1], D[j, j-1] = right*w, -right
        end
    end
    p = try
        filter(isfinite, eigvals(-C, D))
    catch
        _, _, _, _, alpha, beta = schur(complex(-C), complex(D))
        filter(isfinite, alpha ./ beta)
    end
    return RFA.residues(r, p)
end
