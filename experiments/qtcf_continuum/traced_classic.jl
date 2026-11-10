# Record the weights computed by Classic construction without changing its updates.
struct TracedClassic <: RFA.ThieleMethod
    reductions::Vector{NamedTuple}
end
RFA._new_weight(::TracedClassic,z,w,zn,yn) = RFA._new_weight(Classic(),z,w,zn,yn)
function RFA._new_weight(method::TracedClassic, r::QuadraticThiele, z, y)
    h = RFA._new_weight(Classic(), r, z, y)
    g = y
    peak = abs(g)
    for (a, b, pair) in zip(r.apolys, r.bpolys, r.pairs)
        n = RFA._evaluate_block(a, pair, z; numerator=true)
        d = g - RFA._evaluate_block(b, pair, z)
        g = iszero(d) ? oftype(d, iszero(n) ? NaN : Inf) : n/d
        peak = max(peak, abs(g))
    end
    push!(method.reductions, (;block=length(r.pairs)+1, reduced_max=abs(h),
        intermediate_max=peak, slope=NaN, factor=1.0, outcome=isfinite(h) ? "original" : "failed",
        detail=isfinite(h) ? "" : "nonfinite original reduced value"))
    return h
end
