# Independent spectral extraction from the exact reduced stored rational.
# The including driver provides exact_model, peval, and pderiv.

function companion_seeds(q)
    first_nonzero = findfirst(!iszero, q)
    zero_count = first_nonzero - 1
    zero_count <= 1 || throw(ArgumentError("multiple pole at zero"))
    n = length(q) - first_nonzero
    n == 0 && return zeros(eltype(q), zero_count)

    # z = 2^shift*w balances the first and last nonzero coefficients.
    c = q[first_nonzero:end]
    shift = round(Int, (log2(abs(first(c))) - log2(abs(last(c)))) / n)
    unit = ldexp(one(real(last(c))), shift)
    monic = [c[j]/last(c) * ldexp(one(unit), shift*(j-n-1)) for j in 1:n]
    pencil = zeros(eltype(q), n, n)
    for j in 1:n-1
        pencil[j,j+1] = 1
    end
    pencil[end,:] = -monic
    return vcat(zeros(eltype(q), zero_count), unit .* eigvals(pencil))
end

function polish_roots!(z, q, dq)
    tolerance = sqrt(eps(real(one(eltype(q)))))
    for i in eachindex(z)
        for _ in 1:12
            derivative = peval(dq,z[i])
            iszero(derivative) && break
            step = peval(q,z[i])/derivative
            z[i] -= step
            abs(step) <= tolerance*max(one(abs(step)),abs(z[i])) && break
        end
    end
    return z
end

function spectral_at_precision(p_exact, q_exact, precision_bits)
    setprecision(BigFloat, precision_bits) do
        p, q = Complex{BigFloat}.(p_exact), Complex{BigFloat}.(q_exact)
        dq = pderiv(q)
        poles = polish_roots!(companion_seeds(q),q,dq)
        length(poles)==length(q)-1 || error("incomplete pole list")
        all(isfinite,poles) || error("nonfinite candidate pole")
        sort!(poles;by=z -> (real(z),imag(z)))
        bound = abs.(q)
        steps = [abs(peval(q,z)/peval(dq,z)) for z in poles]
        for i in eachindex(poles)
            z = poles[i]
            scale = max(one(abs(z)),abs(z))
            separation = minimum((abs(z-w) for (j,w) in enumerate(poles) if i!=j);init=BigFloat(Inf))
            relative_residual = abs(peval(q,z))/peval(bound,abs(z))
            # A root at zero can give 0/0 in this normalization.
            residual_ok = iszero(z) && iszero(first(q)) || relative_residual < big"1e-45"
            residual_ok && steps[i] < min(big"1e-40"*scale,separation/64) ||
                error("pole iteration did not resolve a simple root")
        end
        residues = [peval(p,z)/peval(dq,z) for z in poles]
        all(isfinite,residues) || error("nonfinite residue")
        return (;poles,residues,precision_bits)
    end
end

# Recompute at higher precision until both poles and residues stabilize.
function repaired_poles(r, kind=:QTCF)
    m = exact_model(r,kind)
    previous = nothing
    for bits in (256,384,640,1024,1536)
        result = try
            spectral_at_precision(m.p,m.q,bits)
        catch e
            e isa Union{ErrorException,LinearAlgebra.LAPACKException} || rethrow()
            continue
        end
        if previous !== nothing
            # Match one to one by location, since conjugate roots need not
            # sort in the same order at different working precisions.
            unused = collect(eachindex(previous.poles))
            consistent = true
            for (z,a) in zip(result.poles,result.residues)
                _, k = findmin(j -> abs(previous.poles[j]-z),unused)
                j = unused[k]
                w,b = previous.poles[j],previous.residues[j]
                consistent &= abs(w-z) <= big"1e-26"*max(1,abs(z))
                consistent &= abs(b-a) <= big"1e-23"*max(big"1e-60",abs(a))
                deleteat!(unused,k)
            end
            consistent && return result
        end
        previous = result
    end
    error("pole and residue extraction did not stabilize")
end
