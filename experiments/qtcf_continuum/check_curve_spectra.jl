# Rebuild curve poles from the exact stored denominator.
# Compare stabilized high-precision residues with the native pencil output.
include("certify_intervals.jl")
include("spectral_reference.jl")
include("check_pole_quality.jl")

upper_abs(z) = abs(real(z))+abs(imag(z))
lower_abs(z) = max(abs(real(z)),abs(imag(z)))

function certify_simple_roots(q, roots)
    dq = pderiv(q)
    ddq = pderiv(dq)
    bound_second = upper_abs.(ddq)
    bound_third = upper_abs.(pderiv(ddq))
    disks = Tuple{EQ,QQ}[]
    for z in roots
        magnitude = max(abs(real(z)),abs(imag(z)))
        exponent = iszero(magnitude) ? 0 : max(0,Base.exponent(magnitude))
        unit = ldexp(big"1.0",exponent-180)
        c = exact(complex(unit*round(BigInt,real(z)/unit),
            unit*round(BigInt,imag(z)/unit)))
        constant = upper_abs(peval(q,c))
        linear = lower_abs(peval(dq,c))
        scale = max(one(QQ),upper_abs(c))
        radius = nothing
        # A dominant linear term certifies one root inside the disk.
        for bits in (55,75,95,115,135,155)
            candidate = scale / (big(1)<<bits)
            bound = peval(bound_second,upper_abs(c)+candidate)*candidate^2/2
            if linear*candidate > constant+bound
                radius = candidate
                break
            end
        end
        if radius === nothing
            second = upper_abs(peval(ddq,c))
            for bits in (55,75,95,115,135,155)
                candidate = scale / (big(1)<<bits)
                bound = second*candidate^2/2 +
                    peval(bound_third,upper_abs(c)+candidate)*candidate^3/6
                if linear*candidate > constant+bound
                    radius = candidate
                    break
                end
            end
        end
        radius === nothing && error("No isolating disk around $(ComplexF64(c))")
        @assert all(abs2(c-w)>(radius+r)^2 for (w,r) in disks)
        push!(disks,(c,radius))
    end
    @assert length(disks)==degree(q)
    return disks
end

function audit_curve_spectra()
    BLAS.set_num_threads(1)
    fits = deserialize(joinpath(DEST,"results.jls"))
    cases = Dict(case.slug=>case for case in make_cases())
    poles = NamedTuple[]
    summaries = NamedTuple[]
    for ((slug,method),fit) in sort!(collect(fits);by=first)
        case = cases[slug]
        (case.family === :interval || fit.best == 0) && continue
        model = exact_model(fit.r, fit.kind)
        origin_pole = iszero(first(model.q))
        rotation = Complex{BigFloat}(case.rotation)
        reflection = z -> rotation*case.reflection(z/rotation)
        diameter = case_metrics(case).diameter*abs(case.rotation)
        exact_symmetry = denominator_has_exact_symmetry(model.q,case.reflection)
        if degree(model.q)==0
            push!(summaries,(;case=slug,method,degree=0,certified_simple_poles=0,
                native_count=0,unique_native_matches=0,native_poles_outside_1e8=0,
                native_residues_outside_1e8=0,native_max_relative_pole_difference=NaN,
                native_max_relative_residue_difference=NaN,
                origin_pole,nearest_pole_real=NaN,nearest_pole_imag=NaN,nearest_residue_real=NaN,
                nearest_residue_imag=NaN,precision_bits=0,
                exact_denominator_reflection_symmetry=exact_symmetry,
                pole_reflection_defect=NaN,pole_reflection_defect_scaled=NaN,
                pole_reflection_defect_error_bound=NaN))
            continue
        end
        reference = repaired_poles(fit.r, fit.kind)
        disks = certify_simple_roots(model.q,reference.poles)
        zp = reference.poles .* rotation
        defect, defect_error = if exact_symmetry
            0.0, 0.0
        else
            certified_centers = [Complex{BigFloat}(c)*rotation for (c,_) in disks]
            # Moving each pole within its certified disk changes each nearest distance
            # by at most twice the largest radius, including the reflected point.
            bound = 2abs(rotation)*maximum(BigFloat(r) for (_,r) in disks)
            Float64(max_symmetry_defect(certified_centers,reflection)), Float64(bound)
        end
        A = reference.residues .* BigFloat(abs2(case.rotation))
        native_z,native_residues = extracted(fit.r, fit.kind)
        native_z .*= case.rotation
        native_residues .*= abs2(case.rotation)
        matches = Int[]
        pole_differences,residue_differences = Float64[],Float64[]
        for (z,a) in zip(native_z,native_residues)
            isfinite(z) || continue
            delta,i = findmin(w->abs(z-w),zp)
            push!(matches,i)
            push!(pole_differences,Float64(delta/max(1,abs(zp[i]))))
            push!(residue_differences,Float64(abs(a-A[i])/abs(A[i])))
        end
        _,nearest = findmin(abs,zp)
        push!(summaries,(;case=slug,method,degree=degree(model.q),certified_simple_poles=length(disks),
            native_count=length(native_z),unique_native_matches=length(unique(matches)),
            native_poles_outside_1e8=count(>(1e-8),pole_differences),
            native_residues_outside_1e8=count(x->!isfinite(x)||x>1e-8,residue_differences),
            native_max_relative_pole_difference=maximum(pole_differences;init=0.0),
            native_max_relative_residue_difference=maximum(residue_differences;init=0.0),
            origin_pole,nearest_pole_real=Float64(real(zp[nearest])),nearest_pole_imag=Float64(imag(zp[nearest])),
            nearest_residue_real=Float64(real(A[nearest])),nearest_residue_imag=Float64(imag(A[nearest])),
            precision_bits=reference.precision_bits,
            exact_denominator_reflection_symmetry=exact_symmetry,
            pole_reflection_defect=defect,pole_reflection_defect_scaled=defect/diameter,
            pole_reflection_defect_error_bound=defect_error))
        for (i,(z,a)) in enumerate(zip(zp,A))
            push!(poles,(;case=slug,method,index=i,real=Float64(real(z)),imag=Float64(imag(z)),
                residue_real=Float64(real(a)),residue_imag=Float64(imag(a)),residue_abs=Float64(abs(a)),
                finite_pole=isfinite(ComplexF64(z)),finite_residue=isfinite(ComplexF64(a)),
                lifted_central_multiple_pole=false,simple_residue_does_not_describe_full_principal_part=false))
        end
        table("verified_curve_poles",poles)
        table("curve_spectral_verification",summaries)
        println(last(summaries))
        flush(stdout)
    end
    expected = count(((key, fit),) -> cases[first(key)].family !== :interval && fit.best > 0, fits)
    @assert length(summaries) == expected
    table("curve_spectral_verification",summaries)
    println("COMPLETED CURVE SPECTRA ",length(summaries)," models; ",length(poles)," simple poles")
end

abspath(PROGRAM_FILE)==(@__FILE__) && audit_curve_spectra()
