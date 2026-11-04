using RationalFunctionApproximation, ComplexRegions, CairoMakie, LaTeXStrings
using Logging, Printf, LinearAlgebra

local_source = joinpath(@__DIR__, "..", "src", "RationalFunctionApproximation.jl")
realpath(pathof(RFA)) == realpath(local_source) || error(
    "These experiments require the local RFA package. Run this file with --project=test from this checkout.")

const EXPERIMENT_MODE = isdefined(@__MODULE__, :QTCF_EXPERIMENT_MODE) ?
    QTCF_EXPERIMENT_MODE : :continuum
EXPERIMENT_MODE in (:continuum, :discrete) || error("Unknown experiment mode: $EXPERIMENT_MODE")

include("qtcf_parameters.jl")

# Interval data and independent validation points.
fine = T(2) .^ ABS_LOG_EXPONENTS
abs_check = sort!(unique(vcat([-one(T) + T(2)*k/(ABS_UNIFORM_N-1)
        for k in 0:ABS_UNIFORM_N-1],
    fine, -fine, fine .- one(T))))
interval_cases = (
    (slug="real", name="|x| on [-1,1]", title=raw"$|x|$ on $[-1,1]$",
        f=x -> abs(real(x)), domain=Segment{T}(-1, 1), points=abs_check),
    (slug="imaginary", name="|z| on [-i,i]", title=raw"$|z|$ on $[-i,i]$",
        f=abs, domain=Segment(-1im, 1im), points=im .* abs_check),
)

# Curves and their analytic derivatives.
function paper_closed_curve(args...)
    curve = Curve(args...)
    point(curve, 0) == point(curve, 1) || error("Curve endpoints must agree exactly")
    return ComplexRegions._closedcurve(curve)
end

paper_ellipse(a, b) = paper_closed_curve(
    t -> a*cospi(2t) + b*im*sinpi(2t),
    t -> 2pi*(-a*sinpi(2t) + b*im*cospi(2t)))

paper_upper_ellipse = Curve(t -> cospi(1-t) + 0.5im*sinpi(1-t),
    t -> pi*sinpi(1-t) - 0.5im*pi*cospi(1-t), 0, 1)
paper_lobed_full = paper_closed_curve(
    t -> (1 + 0.2sinpi(5*(t/pi)))*cispi(t/pi),
    t -> cispi(t/pi)*(cospi(5*(t/pi)) + im*(1 + 0.2sinpi(5*(t/pi)))), -pi/2, 3pi/2)
paper_smooth_s_curve = Curve(t -> -0.75sinpi(2t) + im*(2-4t),
    t -> -1.5pi*cospi(2t) - 4im, 0, 1)
paper_semicircle_s_curve = Curve(
    t -> t <= 0.5 ? im + cispi(0.5+2t) : -im + cispi(0.5-(2t-1)),
    t -> t <= 0.5 ? 2pi*im*cispi(0.5+2t) : -2pi*im*cispi(0.5-(2t-1)), 0, 1)

function super6_point(t)
    c, s = cospi(2t), sinpi(2t)
    cispi(2t) / (c^6 + s^6)^(1/6)
end

function super6_deriv(t)
    c, s = cospi(2t), sinpi(2t)
    g = c^6 + s^6
    2pi*cispi(2t)*(g^(-7/6)*(s*c^5-c*s^5) + im*g^(-1/6))
end

paper_superellipse6 = paper_closed_curve(super6_point, super6_deriv, 0, 1)
paper_inlet = paper_closed_curve(
    t -> INLET_OFFSET + sum(INLET_COS[k]*cospi(k*t/pi) + im*INLET_SIN[k]*sinpi(k*t/pi)
        for k in eachindex(INLET_COS)),
    t -> sum(-k*INLET_COS[k]*sinpi(k*t/pi) + im*k*INLET_SIN[k]*cospi(k*t/pi)
        for k in eachindex(INLET_COS)), 0, 2pi)

qtcf_cases = map((
    ("ellipse_1_025", "ellipse (1, 1/4) Schwarz", paper_ellipse(1, 0.25), conj,
        raw"$z(t)=\cos t+\frac{i}{4}\sin t$", 1.0+0.0im),
    ("ellipse_1_05", "ellipse (1, 1/2) Schwarz", paper_ellipse(1, 0.5), conj,
        raw"$z(t)=\cos t+\frac{i}{2}\sin t$", 1.0+0.0im),
    ("upper_ellipse_1_05", "upper ellipse (1, 1/2) Schwarz", paper_upper_ellipse, conj,
        raw"$z(t)=\cos(\pi(1-t))+\frac{i}{2}\sin(\pi(1-t))$", 1im),
    ("five_lobed", "five-lobed Schwarz", paper_lobed_full, conj,
        raw"$z(t)=(1+0.2\sin 5t)e^{it}$", 1im),
    ("smooth_s", "smooth letter S Schwarz", paper_smooth_s_curve, -,
        raw"$z(t)=-0.75\sin(2\pi t)+i(2-4t)$", 1.0+0.0im),
    ("two_semicircle_s", "two-semicircle S Schwarz", paper_semicircle_s_curve, -,
        raw"$z(t)=i+e^{i(\pi/2+2\pi t)}\ (t\leq1/2),\quad z(t)=-i+e^{i(3\pi/2-2\pi t)}\ (t>1/2)$", 1.0+0.0im),
    ("superellipse_6", "superellipse x^6+y^6=1 Schwarz", paper_superellipse6, conj,
        raw"$z(t)=e^{2\pi it}/(\cos^6(2\pi t)+\sin^6(2\pi t))^{1/6}$", 1.0+0.0im),
    ("analytic_inlet", "analytic right inlet Schwarz", paper_inlet, conj,
        raw"$z(t)=0.1-0.1\cos t-0.7\cos2t-0.1\cos3t-i(0.8\sin t+0.4\sin2t+0.2\sin3t)$", 1.0+0.0im),
)) do (slug, name, domain, reflection, formula, rotation)
    constructor = RFA.isclosed(domain) ? paper_closed_curve : Curve
    fit_domain = rotation == 1 ? domain : constructor(
        t -> point(domain, t)/rotation,
        t -> ComplexRegions.tangent(domain, t)/rotation, 0, 1)
    (; slug, name, domain, fit_domain, reflection, function_reflection=reflection, formula, rotation)
end

include("qtcf_output.jl")

# The scalar TCF interface reads its recurrence selectors at fit time.
set_eval_method(EXPERIMENT_EVAL_METHOD)
set_weight_method(EXPERIMENT_WEIGHT_METHODS.TCF)

function fit_experiment(f, domain, method; kwargs...)
    if method === :QTCF
        return approximate(f, domain, QTCF(); weight_method=EXPERIMENT_WEIGHT_METHODS.QTCF,
            eval_method=EXPERIMENT_EVAL_METHOD, kwargs...)
    end
    return approximate(f, domain, method === :AAA ? AAA() : TCF(); kwargs...)
end

# Use the same parameter-midpoint startup in both sampling modes.
qtcf_startup(domain) = (; initial_node=point(domain, length(domain)*QTCF_START_PARAMETER))

# Check one scalar startup followed by exact symmetry pairs.
function check_qtcf_blocks(result, domain, reflection)
    fixed = qtcf_startup(domain).initial_node
    @assert reflection(fixed) == fixed
    for (k, record) in enumerate(result.history)
        r = record.interpolant
        pairs = r.pairs
        @assert length(pairs) == k && length(nodes(r)) == 2k - 1
        @assert first(pairs) == (fixed, fixed)
        for (z, w) in @view pairs[2:end]
            @assert z != w && reflection(z) == w
        end
        @assert all(b -> all(isfinite, b), RFA.weights(r))
    end
end

# Fixed discrete grids share all continuum supports and reflected samples.
continuum_nodes(result) =
    unique([z for record in result.history for z in nodes(record.interpolant)])

# Signed zeros denote the same sample coordinate.
sample_coordinate(x::Real) = iszero(x) ? zero(x) : x
sample_coordinate(z::Complex) = complex(sample_coordinate(real(z)), sample_coordinate(imag(z)))

function support_parameters(domain, result, support)
    mesh_s, mesh_z = collect(result.path, :nodes)
    positions = Dict{eltype(mesh_z),eltype(mesh_s)}()
    for (s, z) in zip(mesh_s, mesh_z)
        get!(positions, sample_coordinate(z), s)
    end
    map(support) do z
        # Center neighbors at the recorded parameter, or the closest mesh point.
        s = get(positions, sample_coordinate(z)) do
            domain isa Segment && return T(ComplexRegions.arg(domain, z))
            mesh_s[argmin(k -> abs(mesh_z[k] - z), eachindex(mesh_z))]
        end
        T(s) / length(domain)
    end
end

function discrete_grid_samples(domain, seeds)
    closed = RFA.isclosed(domain)
    clean(s) = closed ? mod(s, one(T)) : clamp(s, zero(T), one(T))
    mapped(s) = point(domain, length(domain) * clean(s))
    Z = typeof(mapped(zero(T)))
    supports = map(continuum_nodes, seeds)
    support = unique(Z[z for points in supports for z in points])
    parameters = unique(vcat(support_parameters.(Ref(domain), seeds, supports)...))
    # Uniform samples, continuum nodes, and logarithmic parameter neighbors.
    uniform = mapped.(range(zero(T); step=one(T)/DISCRETE_UNIFORM_N,
        length=DISCRETE_UNIFORM_N + !closed))
    neighbors = Z[mapped(s + sign * 10.0^exponent)
        for s in parameters for exponent in DISCRETE_LOG_EXPONENTS for sign in (-1, 1)]
    (; support, parameters, uniform, neighbors, seed_node_counts=length.(supports))
end

function shared_discrete_grid(domain, reflection, seeds;
    samples=discrete_grid_samples(domain, seeds),
    source="uniform+all_recorded_continuum_support_union+parameter_log_neighbors+exact_reflection_closure")
    (; support, parameters, uniform, neighbors, seed_node_counts) = samples
    points = unique(sample_coordinate.(vcat(support, uniform, neighbors)))
    Z = eltype(points)
    base_count = length(points)
    unique!(append!(points, sample_coordinate.(Z.(reflection.(points)))))
    exact = Set(points)
    reflection_closed = all(z -> sample_coordinate(Z(reflection(z))) in exact, points)
    reflection_closed || error("The shared fitting grid is not closed under reflection")
    println("  shared fixed grid: ", base_count, " base points, ", length(points),
        " after exact reflection closure (", length(unique(sample_coordinate.(support))), " continuum interpolation nodes)")
    flush(stdout)
    (; points, exact, parameters, closed=RFA.isclosed(domain), source, base_count,
        reflection_closed, uniform_count=length(unique(sample_coordinate.(uniform))),
        support_count=length(unique(sample_coordinate.(support))), support_parameter_count=length(parameters),
        neighbor_count=length(unique(sample_coordinate.(neighbors))), seed_node_counts)
end

# Preserve the original TCF fitting samples and their union order.
function paper_interval_grid(case)
    logarithmic = T(2) .^ ABS_FITTING_LOG_EXPONENTS
    half = (ABS_UNIFORM_N - 1) ÷ 2
    uniform = collect(-half:half) / T(half)
    realpoints = union(collect(logarithmic), collect(-logarithmic),
        collect(uniform), collect(logarithmic .- one(T)))
    points = case.slug == "real" ? realpoints : im .* realpoints
    reflection = symmetric(case.domain).reflection
    exact = Set(sample_coordinate.(points))
    (; points, exact, closed=false,
        source="original_tcf_ordered_union+uniform_10001+log2_exponents_minus52_to_minus1+left_endpoint_cluster",
        base_count=length(points),
        reflection_closed=all(z -> sample_coordinate(reflection(z)) in exact, points),
        uniform_count=length(uniform), uniform_subdivisions=ABS_UNIFORM_N-1,
        support_count=0, support_parameter_count=0, neighbor_count=0,
        seed_node_counts=(0, 0))
end

function run_intervals(metadata)
    map(interval_cases) do case
        reflection = symmetric(case.domain).reflection
        if EXPERIMENT_MODE === :discrete
            println("\nDiscrete interval case: ", case.name)
            grid = paper_interval_grid(case)
            tcf = fit_experiment(case.f, grid.points, :TCF;
                experiment_parameters(:interval, :TCF; discrete=true)...)
            qtcf = fit_experiment(symmetric(case.f), symmetric(grid.points; reflection), :QTCF;
                experiment_parameters(:interval, :QTCF; discrete=true)...,
                qtcf_startup(case.domain)...)
            for (method, result) in zip((:TCF, :QTCF), (tcf, qtcf))
                record_discrete_run!(metadata, case.slug, method, result, grid, 0,
                    length(validation_grid(case).points))
            end
        else
            println("\nContinuum interval case: ", case.name)
            tcf = fit_experiment(case.f, case.domain, :TCF;
                experiment_parameters(:interval, :TCF)...)
            qtcf = fit_experiment(symmetric(case.f), symmetric(case.domain), :QTCF;
                experiment_parameters(:interval, :QTCF)...,
                parameter_reflection=qtcf_parameter_reflection,
                qtcf_startup(case.domain)...)
        end
        check_qtcf_blocks(qtcf, case.domain, reflection)
        (; case..., runs=(tcf, qtcf))
    end
end

function curve_points(domain, n; periodic=false)
    parameters = range(0, length(domain); length=n + periodic)
    # Closed-curve validation counts the common endpoint only once.
    periodic && (parameters = parameters[1:end-1])
    return [point(domain, t) for t in parameters]
end

function fit_schwarz(case, method; points=nothing)
    domain = points === nothing ? case.fit_domain : points
    startup = if method === :QTCF
        options = qtcf_startup(case.fit_domain)
        points === nothing ? (; options...,
            parameter_reflection=qtcf_parameter_reflection) : options
    else
        (;)
    end
    f, domain = method === :QTCF ?
        (symmetric(conj; reflection=case.function_reflection),
            symmetric(domain; reflection=case.reflection)) : (conj, domain)
    result = fit_experiment(f, domain, method;
        experiment_parameters(:schwarz, method; discrete=points !== nothing)..., startup...)
    if method === :QTCF
        check_qtcf_blocks(result, case.fit_domain, case.reflection)
    end
    result
end

function run_method(case, method, audit, target, metadata; grid=nothing, seed_nodes=0)
    start = time_ns()
    result = fit_schwarz(case, method; points=grid === nothing ? nothing : grid.points)
    elapsed = (time_ns() - start) / 1e9
    rational = get_function(result)
    extraction_start = time_ns()
    pole, residue = residues(rational)
    extraction_seconds = (time_ns() - extraction_start) / 1e9
    length(pole) == length(residue) || error("Pole and residue counts must agree")
    finite_poles = count(isfinite, pole)
    nonfinite_poles = length(pole) - finite_poles
    nonfinite_residues = count(!isfinite, residue)
    extraction_status = iszero(nonfinite_poles + nonfinite_residues) ? :finite : :nonfinite
    rotation = case.rotation
    pole .*= rotation
    residue .*= abs2(rotation)
    error = abs(rotation) * max_error(rational, target, audit.points)
    rational_degrees = degrees(rational)
    stop = status(result)
    grid === nothing || record_discrete_run!(metadata, case.slug, method, result, grid,
        seed_nodes, length(audit.points))
    (; pole, residue, error, elapsed, extraction_seconds, extraction_status,
        finite_poles, nonfinite_poles, nonfinite_residues, rational_degrees, failed=false,
        selected_nodes=length(nodes(rational)), recorded=length(result.history),
        returned_record=stop.best, sampled_error=stop.error, stop_reason=stop.reason,
        validation_points=length(audit.points), validation_nodes=0,
        validation_grid_sha256=audit.metadata.grid_sha256)
end

function run_case(case, metadata)
    println("\nSchwarz case: ", case.name)
    grid = if EXPERIMENT_MODE === :discrete
        seeds = map(METHODS) do method
            println("  continuum grid seed ", method)
            fit_schwarz(case, method)
        end
        shared_discrete_grid(case.fit_domain, case.reflection, seeds)
    end
    audit = validation_grid(case)
    target = conj.(audit.points)
    rows = Dict{Symbol,Any}()
    for (j, method) in enumerate(METHODS)
        print("  ", EXPERIMENT_MODE, " ", method, " ... ")
        flush(stdout)
        try
            row = run_method(case, method, audit, target, metadata; grid,
                seed_nodes=grid === nothing ? 0 : grid.seed_node_counts[j])
            rows[method] = row
            @printf("error %.3e, degree %s, %.3f s\n",
                row.error, string(row.rational_degrees), row.elapsed)
        catch exception
            exception isa InterruptException && rethrow()
            rows[method] = (
                pole=ComplexF64[], residue=ComplexF64[],
                error=NaN, elapsed=NaN, extraction_seconds=NaN, extraction_status=:failed,
                finite_poles=0, nonfinite_poles=0, nonfinite_residues=0,
                rational_degrees=(-1, -1), failed=true,
                selected_nodes=0, recorded=0, returned_record=0,
                sampled_error=NaN, stop_reason=:failed,
                validation_points=0, validation_nodes=0, validation_grid_sha256="")
            println("failed: ", exception)
        end
        flush(stdout)
    end
    (; case, curve=curve_points(case.domain, SCHWARZ_CURVE_N), rows)
end

function run_experiments()
    capture_run_provenance!()
    global_logger(ConsoleLogger(stderr, Logging.Error))
    BLAS.set_num_threads(1)
    println("RationalFunctionApproximation ", pkgversion(RFA))
    println("Package source: ", pathof(RFA))
    println("Running all paper experiments.")
    EXPERIMENT_MODE === :discrete && println(
        "Discrete mode: original TCF interval samples and continuum-supported Schwarz grids.")
    println("Output: ", OUTPUT)

    metadata = NamedTuple[]
    interval_results = run_intervals(metadata)
    schwarz_results = [run_case(case, metadata) for case in qtcf_cases]
    write_outputs(interval_results, schwarz_results, metadata)
    failed_runs = [(result.case.slug, method) for result in schwarz_results
        for method in METHODS if result.rows[method].failed]
    isempty(failed_runs) || error("Some Schwarz experiments failed; inspect numbers/schwarz_summary.csv: $failed_runs")
    println("Completed all 28 fits.")
    println("Saved figures and numerical data to ", OUTPUT)
    (; interval_results, schwarz_results, metadata)
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiments()
