# Certify reduced denominator degrees from results.jls without refitting.
# Pass --histories to include every recorded prefix.
include("certify_intervals.jl")

# p is prime and -1 is nonsquare, so Gaussian dyadic coefficients reduce
# to F_{p²}. The reduction is exact because p does not divide a power of 2.
const MODULUS = 65519
const FF = Complex{Int}
const DYADIC_POWERS = [powermod(e >= 0 ? 2 : invmod(2, MODULUS), abs(e), MODULUS)
    for e in -1200:1200]

red(x::FF) = complex(mod(real(x), MODULUS), mod(imag(x), MODULUS))
function real_residue(x::Real)
    isfinite(x) || error("Cannot certify a nonfinite coefficient")
    f, exponent = frexp(Float64(x))
    shift = exponent - 53
    mantissa = Int(ldexp(f, 53))
    mod(mod(mantissa, MODULUS) * DYADIC_POWERS[shift + 1201], MODULUS)
end
residue(x) = complex(real_residue(real(x)), real_residue(imag(x)))
function field_inverse(x::FF)
    norm_square = mod(abs2(x), MODULUS)
    red(conj(x) * invmod(norm_square, MODULUS))
end

function modular_product(a, b)
    c = zeros(FF, length(a) + length(b) - 1)
    for j in eachindex(b), i in eachindex(a)
        c[i+j-1] = red(c[i+j-1] + a[i]*b[j])
    end
    trim(c)
end

function modular_sum(a, b)
    c = zeros(FF, max(length(a), length(b)))
    c[1:length(a)] .= a
    for i in eachindex(b)
        c[i] = red(c[i]+b[i])
    end
    trim(c)
end

function modular_blocks(r)
    d = r isa RFA.QuadraticThiele ? getfield(r, :data) : r
    if d isa Thiele
        return [[residue(b)] for b in weights(d)], [[residue(-z), one(FF)] for z in nodes(d)]
    end
    b = [residue.(v) for v in d.bpolys]
    a = Vector{FF}[]
    for (v, pair) in zip(d.apolys, d.pairs)
        first_factor = [residue(-pair[1]), one(FF)]
        product = if length(v) == 2
            first_factor
        else
            modular_product(first_factor, [residue(-pair[2]), one(FF)])
        end
        push!(a, red.(residue(last(v)) .* product))
    end
    b, a
end

function modular_continuants(r, kind)
    if r isa RFA.QuadraticThiele && r.chart !== nothing
        @assert kind === :QTCF
        return chart_continuants(r,residue,modular_sum,modular_product)
    end
    b, a = modular_blocks(r)
    isempty(b) && return FF[0], FF[1]
    prevp, prevq, p, q = FF[1], FF[0], trim(first(b)), FF[1]
    for j in 2:length(b)
        p, q, prevp, prevq =
            modular_sum(modular_product(b[j], p), modular_product(a[j-1], prevp)),
            modular_sum(modular_product(b[j], q), modular_product(a[j-1], prevq)), p, q
    end
    if kind in (:even, :odd)
        p, q = lift_poly(p), lift_poly(q)
        kind === :odd && (p = vcat(FF(0), p))
    end
    trim(p), trim(q)
end

function modular_remainder(a, b)
    r = copy(a)
    inv_leading = field_inverse(last(b))
    for shift in (length(r)-length(b)):-1:0
        lead = red(r[end] * inv_leading)
        for j in eachindex(b)
            r[j+shift] = red(r[j+shift] - lead*b[j])
        end
        pop!(r)
    end
    isempty(r) ? FF[0] : trim(r)
end

function modular_gcd_degree(a, b)
    while !iszero_poly(b)
        a, b = b, modular_remainder(a, b)
    end
    degree(a)
end

modular_value(a, z) = foldr((c, v) -> red(c+z*v), a; init=zero(FF))

function preceding_nodes(r)
    d = r isa RFA.QuadraticThiele ? getfield(r, :data) : r
    d isa Thiele && return nodes(d)[1:end-1]
    roots = ComplexF64[]
    for j in 1:length(d.bpolys)-1
        push!(roots, d.pairs[j][1])
        length(d.apolys[j]) == 3 && push!(roots, d.pairs[j][2])
    end
    roots
end

function exact_pair_value(b, a, z)
    prevp, prevq, p, q = one(EQ), zero(EQ), peval(first(b), z), one(EQ)
    for j in 2:length(b)
        bj, aj = peval(b[j], z), peval(a[j-1], z)
        p, q, prevp, prevq = bj*p+aj*prevp, bj*q+aj*prevq, p, q
    end
    p, q
end

function exact_cancelled_roots(r, p, q)
    roots = preceding_nodes(r)
    allunique(roots) || return nothing
    # The continuant determinant has only these distinct linear factors.
    # Use the finite field to identify possible roots, then evaluate the
    # exact projective pair at them to rule out modular false positives.
    candidates = filter(roots) do z
        modz = residue(z)
        iszero(modular_value(p, modz)) && iszero(modular_value(q, modz))
    end
    b, a = exact_blocks(r)
    count(candidates) do z
        all(iszero, exact_pair_value(b, a, exact(z)))
    end
end

function certify_degree(r, kind; exact_fallback)
    bound = last(original_degrees(r, kind))
    p, q = modular_continuants(r, kind)
    qmod = iszero_poly(q) ? -1 : degree(q)
    gcdmod = modular_gcd_degree(p, q)
    # The known upper bound and full modular degree certify the exact
    # unreduced degree of q. Since q's leading coefficient is a unit at
    # p, a positive-degree exact common factor would also reduce with
    # positive degree. A modular gcd of degree zero therefore proves
    # that q has no cancellations.
    if qmod == bound && gcdmod == 0
        return (; nominal_degree=bound, reduced_degree=bound,
            modular_degree=qmod, modular_gcd_degree=gcdmod, status="certified_modular")
    elseif exact_fallback
        chart = r isa RFA.QuadraticThiele && r.chart !== nothing
        cancelled = if chart && qmod==bound && gcdmod==1 && iszero(first(p)) && iszero(first(q))
            # Both exact constant terms vanish, so z divides the exact pair.
            # The degree-one modular gcd bounds all exact common factors.
            all(iszero,exact_chart_pair(r,zero(EQ))) ? 1 : nothing
        elseif !chart && qmod == bound && kind in (:QTCF, :forced, :native)
            exact_cancelled_roots(r, p, q)
        end
        reduced = if cancelled === nothing
            m = exact_model(r, kind)
            @assert !iszero_poly(m.q)
            degree(m.q)
        else
            bound - cancelled
        end
        @assert 0 <= reduced <= bound
        return (; nominal_degree=bound, reduced_degree=reduced,
            modular_degree=qmod, modular_gcd_degree=gcdmod, status="certified_exact")
    end
    (; nominal_degree=bound, reduced_degree=-1,
        modular_degree=qmod, modular_gcd_degree=gcdmod, status="requires_exact_check")
end

function main(args)
    certify_histories = args == ["--histories"]
    isempty(args) || certify_histories || error("Pass --histories to certify all recorded models")
    @assert all(MODULUS % i != 0 for i in 2:isqrt(MODULUS))
    @assert MODULUS % 4 == 3
    fits = deserialize(joinpath(DEST, "results.jls"))
    rows = NamedTuple[]
    for ((slug, method), fit) in sort!(collect(fits); by=first)
        if certify_histories
            fit.best == 0 && continue
            for (block, record) in enumerate(fit.history)
                row = certify_degree(record.interpolant, fit.kind; exact_fallback=true)
                push!(rows, (;case=slug, method, block, returned=block==fit.best, row...))
            end
            println(slug, " ", method, " completed ", length(fit.history), " records")
        elseif fit.best==0
            push!(rows, (;case=slug, method, nominal_degree=-1, reduced_degree=-1,
                modular_degree=-1, modular_gcd_degree=-1, status="no_returned_fit"))
        else
            row = certify_degree(fit.r, fit.kind; exact_fallback=true)
            push!(rows, (;case=slug, method, row...))
            println(last(rows))
        end
        flush(stdout)
    end
    filename = certify_histories ? "reduced_history_degrees.csv" : "reduced_returned_degrees.csv"
    save_cert(@__DIR__, filename, rows)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
