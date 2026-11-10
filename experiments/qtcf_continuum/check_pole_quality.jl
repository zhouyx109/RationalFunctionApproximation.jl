# Check native pencil candidates against the exact stored denominator.
if !isdefined(@__MODULE__, :exact_model)
    include("certify_intervals.jl")
end
using DelimitedFiles: readdlm

const LARGE_POLE_RESIDUAL = 1e-8

has_simple_small_denominator(q) = degree(q)==1 ||
    (degree(q)==2 && !iszero(q[2]^2-4q[1]*q[3]))

function load_native_poles()
    data, header = readdlm(joinpath(DEST,"poles.csv"), ',', String, '\n'; header=true)
    keys = Tuple(Symbol.(vec(header)))
    groups = Dict{Tuple{String,String},Vector{NamedTuple}}()
    for row in eachrow(data)
        r = NamedTuple{keys}(Tuple(row))
        push!(get!(Vector{NamedTuple},groups,(r.case,r.method)),r)
    end
    return groups
end

# This is a residual of the stored polynomial, not a bound on pole location.
function pole_residual(q, magnitudes, z)
    value = abs(peval(q,z))
    scale = peval(magnitudes,abs(z))
    return iszero(scale) ? (iszero(value) ? 0.0 : Inf) : Float64(value/scale)
end

# Test polynomial reflection up to a nonzero constant factor.
function denominator_has_exact_symmetry(q, reflection)
    if reflection(one(EQ))==one(EQ)
        @assert reflection(EQ(im))==-EQ(im)
        leading = last(q)
        return all(c->c*conj(leading)==conj(c)*leading,q)
    end
    @assert reflection(one(EQ))==-one(EQ) && reflection(EQ(im))==-EQ(im)
    return all(j->iseven(length(q)-j) || iszero(q[j]),eachindex(q))
end

# For the low-degree returns missed by the pencil, solve the exact reduced
# denominator directly. Compute the smaller quadratic root using their product.
function small_denominator_roots(q)
    degree(q)==1 && return [Complex{BigFloat}(-q[1]/q[2])]
    @assert degree(q)==2
    c,b,a = q
    discriminant = b*b-4a*c
    @assert !iszero(discriminant)
    d = sqrt(Complex{BigFloat}(discriminant))
    plus,minus = -b+d,-b-d
    large = (abs(plus)>=abs(minus) ? plus : minus)/(2a)
    return [large,c/(a*large)]
end

function recover_odd_poles(m, case, method)
    pp,dq = Complex{BigFloat}.(m.p),Complex{BigFloat}.(pderiv(m.q))
    roots = small_denominator_roots(m.q)
    rotation = Complex{BigFloat}(case.rotation)
    rows = NamedTuple[]
    for (i,z) in enumerate(roots)
        @assert !iszero(peval(dq,z))
        physical_pole = ComplexF64(rotation*z)
        residue = ComplexF64(abs2(rotation)*peval(pp,z)/peval(dq,z))
        push!(rows,(;case=case.slug,method,index=i,real=real(physical_pole),imag=imag(physical_pole),
            residue_real=real(residue),residue_imag=imag(residue),residue_abs=abs(residue),
            finite_pole=isfinite(physical_pole),finite_residue=isfinite(residue),
            lifted_central_multiple_pole=false,simple_residue_does_not_describe_full_principal_part=false))
    end
    return rows
end

function check_nonfinite_validation(fit, m, case, method)
    bad = [!isfinite(lifted(fit.r,z,fit.kind)) for z in case.check]
    any(bad) || return nothing
    p,q = Complex{BigFloat}.(m.p),Complex{BigFloat}.(m.q)
    poles, peak, peak_on_bad = 0,0.0,0.0
    for (z,stock_nonfinite) in zip(case.check,bad)
        w = Complex{BigFloat}(z)
        d = peval(q,w)
        value = if stock_nonfinite || iszero(d)
            # Resolve stock failures on the exact dyadic stored rational.
            v = exact(z)
            denominator = peval(m.q,v)
            if iszero(denominator)
                poles+=1
                Complex{BigFloat}(Inf)
            else
                Complex{BigFloat}(peval(m.p,v)/denominator)
            end
        else
            peval(p,w)/d
        end
        err = Float64(abs(case.f(z)-value))
        peak = max(peak,err)
        stock_nonfinite && (peak_on_bad = max(peak_on_bad,err))
    end
    return (;case=case.slug,method,precision_bits=precision(BigFloat),points=length(bad),
        stock_nonfinite=count(identity,bad),exact_poles_on_validation_points=poles,
        stored_rational_sampled_error=peak,stored_rational_error_on_stock_nonfinite_points=peak_on_bad)
end

function check_pole_quality()
    fits = deserialize(joinpath(DEST,"results.jls"))
    native = load_native_poles()
    checks,recovered,nonfinite_checks = NamedTuple[],NamedTuple[],NamedTuple[]
    setprecision(BigFloat,256) do
        cases = Dict(case.slug=>case for case in make_cases())
        for ((slug,method),fit) in sort!(collect(fits);by=first)
            case = cases[slug]
            key = (slug,method)
            if fit.best==0
                push!(checks,(;case=slug,method,status="no_returned_fit",exact_reduced_degree=-1,
                    exact_denominator_reflection_symmetry=missing,native_finite_candidates=0,
                    native_finite_residues=0,native_reflection_defect=missing,
                    max_native_denominator_residual=missing,native_residual_above_1e_minus_8=0,
                    recovered_finite_poles=0,recovered_reflection_defect=missing))
                continue
            end
            m = exact_model(fit.r,fit.kind)
            if slug in ("real","imaginary","odd_sine")
                reevaluation = check_nonfinite_validation(fit,m,case,method)
                reevaluation===nothing || push!(nonfinite_checks,reevaluation)
            end
            qq = Complex{BigFloat}.(m.q)
            magnitudes = abs.(qq)
            rows = get(native,key,NamedTuple[])
            candidates = [complex(parse(Float64,r.real),parse(Float64,r.imag))
                for r in rows if r.finite_pole=="true"]
            rotation = Complex{BigFloat}(case.rotation)
            reflection = z -> case.rotation*case.reflection(z/case.rotation)
            errors = [pole_residual(qq,magnitudes,Complex{BigFloat}(z)/rotation) for z in candidates]
            recovered_rows = NamedTuple[]
            if case.parity===:odd && length(candidates)<degree(m.q) && has_simple_small_denominator(m.q)
                recovered_rows = recover_odd_poles(m,case,method)
                append!(recovered,recovered_rows)
            end
            recovered_candidates = [complex(r.real,r.imag) for r in recovered_rows]
            @assert all(z -> pole_residual(qq,magnitudes,Complex{BigFloat}(z)/rotation)<1e-12,
                recovered_candidates)
            push!(checks,(;case=slug,method,status="returned_fit",exact_reduced_degree=degree(m.q),
                exact_denominator_reflection_symmetry=denominator_has_exact_symmetry(m.q,case.reflection),
                native_finite_candidates=length(candidates),
                native_finite_residues=count(r->r.finite_residue=="true",rows),
                native_reflection_defect=max_symmetry_defect(candidates,reflection),
                max_native_denominator_residual=isempty(errors) ? missing : maximum(errors),
                native_residual_above_1e_minus_8=count(>(LARGE_POLE_RESIDUAL),errors),
                recovered_finite_poles=length(recovered_candidates),
                recovered_reflection_defect=max_symmetry_defect(recovered_candidates,reflection)))
            println(last(checks))
            flush(stdout)
        end
    end
    table("pole_quality",checks)
    table("recovered_poles",recovered)
    table("nonfinite_validation_recheck",nonfinite_checks)
end

abspath(PROGRAM_FILE)==(@__FILE__) && check_pole_quality()
