using RationalFunctionApproximation, ComplexRegions, CairoMakie, LaTeXStrings
using Logging, Printf, LinearAlgebra

isdefined(RFA, :QuadraticThiele) || error(
    "These experiments require the local QTCF extension. Run this file with --project=test from this checkout.")

const EXPERIMENT_MODE = isdefined(@__MODULE__, :QTCF_EXPERIMENT_MODE) ?
    QTCF_EXPERIMENT_MODE : :continuum
EXPERIMENT_MODE in (:continuum, :discrete) || error("Unknown experiment mode: $EXPERIMENT_MODE")

include("qtcf_parameters.jl")

# Interval data and independent validation points.
fine = T(2) .^ ABS_LOG_EXPONENTS
abs_check = sort!(unique(vcat(range(-one(T), one(T), ABS_UNIFORM_N),
    fine, -fine, fine .- one(T))))
interval_cases = (
    (slug="real", name="|x| on [-1,1]", title=raw"$|x|$ on $[-1,1]$",
        f=x -> abs(real(x)), domain=Segment{T}(-1, 1), points=abs_check),
    (slug="imaginary", name="|z| on [-i,i]", title=raw"$|z|$ on $[-i,i]$",
        f=abs, domain=Segment(-1im, 1im), points=im .* abs_check),
)

# Curves and their analytic derivatives.
paper_upper_ellipse = Curve(t -> cospi(1-t) + 0.5im*sinpi(1-t),
    t -> pi*sinpi(1-t) - 0.5im*pi*cospi(1-t), 0, 1)
paper_lobed_full = ClosedCurve(
    t -> (1 + 0.2sinpi(5*(t/pi)))*cispi(t/pi),
    t -> cispi(t/pi)*(cospi(5*(t/pi)) + im*(1 + 0.2sinpi(5*(t/pi)))), 0, 2pi)
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

paper_superellipse6 = ClosedCurve(super6_point, super6_deriv, 0, 1)
paper_inlet = ClosedCurve(
    t -> INLET_OFFSET + sum(INLET_COS[k]*cospi(k*t/pi) + im*INLET_SIN[k]*sinpi(k*t/pi)
        for k in eachindex(INLET_COS)),
    t -> sum(-k*INLET_COS[k]*sinpi(k*t/pi) + im*k*INLET_SIN[k]*cospi(k*t/pi)
        for k in eachindex(INLET_COS)), 0, 2pi)

qtcf_cases = map((
    ("ellipse_1_025", "ellipse (1, 1/4) Schwarz", Shapes.ellipse(1, 0.25), conj,
        raw"$z(t)=\cos t+\frac{i}{4}\sin t$", 1.0+0.0im),
    ("ellipse_1_05", "ellipse (1, 1/2) Schwarz", Shapes.ellipse(1, 0.5), conj,
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
    constructor = RFA.isclosed(domain) ? ClosedCurve : Curve
    fit_domain = rotation == 1 ? domain : constructor(
        t -> point(domain, t)/rotation,
        t -> ComplexRegions.tangent(domain, t)/rotation, 0, 1)
    (; slug, name, domain, fit_domain, reflection, function_reflection=reflection, formula, rotation)
end

include("qtcf_output.jl")

# Recurrence selectors must be defined before the numerical calls are compiled.
set_eval_method(OneDiv())
set_weight_method(OneDiv())

# Fixed discrete grids share all continuum supports and reflected samples.
function support_parameters(domain, result)
    mesh_s, mesh_z = collect(result.path, :nodes)
    map(nodes(get_function(result))) do z
        distance, k = findmin(abs.(mesh_z .- z))
        s = T(mesh_s[k])
        # Near a nonsmooth parametrization, parameter roundoff can exceed _same.
        h = 4eps(T) * max(one(T), abs(s))
        bound = maximum(abs(point(domain, RFA._clean_param(domain, s + δ)) - mesh_z[k])
            for δ in (-h, h))
        if !RFA._same(mesh_z[k], z) && distance > bound
            s = RFA._reflection_param(domain, zero(T), _ -> z)
        end
        s / length(domain)
    end
end

function discrete_grid_samples(domain, seeds)
    closed = RFA.isclosed(domain)
    clean(s) = closed ? mod(s, one(T)) : clamp(s, zero(T), one(T))
    mapped(s) = point(domain, length(domain) * clean(s))
    Z = typeof(mapped(zero(T)))
    support = Z[z for result in seeds for z in nodes(get_function(result))]
    parameters = unique(vcat(support_parameters.(Ref(domain), seeds)...))
    # Match the stable-evaluation grid: uniform points, nodes, and log neighbors.
    uniform = mapped.(range(zero(T); step=one(T)/DISCRETE_UNIFORM_N,
        length=DISCRETE_UNIFORM_N))
    neighbors = Z[mapped(s + sign * 10.0^exponent)
        for s in parameters for exponent in DISCRETE_LOG_EXPONENTS for sign in (-1, 1)]
    (; support, parameters, uniform, neighbors)
end

function shared_discrete_grid(domain, reflection, seeds;
    samples=discrete_grid_samples(domain, seeds),
    source="uniform+continuum_support_union+parameter_log_neighbors+exact_reflection_closure")
    (; support, parameters, uniform, neighbors) = samples
    # Signed zeros denote the same coordinate, not two distinct samples.
    cleanzero(x) = iszero(x) ? zero(x) : x
    coordinate(z) = z isa Real ? cleanzero(z) :
        complex(cleanzero(real(z)), cleanzero(imag(z)))
    points = unique(coordinate.(vcat(support, uniform, neighbors)))
    Z = eltype(points)
    base_count = length(points)
    unique!(append!(points, coordinate.(Z.(reflection.(points)))))
    exact = Set(points)
    reflection_closed = all(z -> coordinate(Z(reflection(z))) in exact, points)
    reflection_closed || error("The shared fitting grid is not closed under reflection")
    println("  shared fixed grid: ", base_count, " base points, ", length(points),
        " after exact reflection closure (", length(unique(coordinate.(support))), " continuum support points)")
    flush(stdout)
    (; points, exact, parameters, closed=RFA.isclosed(domain), source, base_count,
        reflection_closed, uniform_count=length(unique(coordinate.(uniform))),
        support_count=length(unique(coordinate.(support))), support_parameter_count=length(parameters),
        neighbor_count=length(unique(coordinate.(neighbors))),
        seed_node_counts=length.(nodes.(get_function.(seeds))))
end

function run_intervals(metadata)
    interval_seeds = map(interval_cases) do case
        println("\nInterval continuum seed: ", case.name)
        tcf = approximate(case.f, case.domain, TCF(); experiment_parameters(:interval, :TCF)...)
        qtcf = approximate(symmetric(case.f), symmetric(case.domain);
            experiment_parameters(:interval, :QTCF)...)
        (tcf, qtcf)
    end

    interval_grids = if EXPERIMENT_MODE === :discrete
        samples = map((case, seeds) -> discrete_grid_samples(case.domain, seeds),
            interval_cases, interval_seeds)
        # Use one coordinate grid for both intervals, so rotation changes no samples.
        common = (; (name => unique(vcat(real.(getproperty(samples[1], name)),
            imag.(getproperty(samples[2], name)))) for name in (:support, :uniform, :neighbors))...)
        parameters = unique(vcat(getproperty.(samples, :parameters)...))
        map(eachindex(interval_cases)) do j
            case, seeds = interval_cases[j], interval_seeds[j]
            rotate = j == 1 ? identity : x -> im * x
            mapped = (; (name => rotate.(getproperty(common, name))
                for name in (:support, :uniform, :neighbors))..., parameters)
            shared_discrete_grid(case.domain, symmetric(case.domain).reflection, seeds;
                samples=mapped,
                source="uniform+all_interval_continuum_support_union+parameter_log_neighbors+exact_reflection_closure")
        end
    end

    map(eachindex(interval_cases)) do j
        case = interval_cases[j]
        tcf, qtcf = interval_seeds[j]
        if EXPERIMENT_MODE === :discrete
            grid = interval_grids[j]
            reflection = symmetric(case.domain).reflection
            println("\nDiscrete interval case: ", case.name)
            tcf = approximate(case.f, grid.points, TCF();
                experiment_parameters(:interval, :TCF; discrete=true)...)
            qtcf = approximate(symmetric(case.f), symmetric(grid.points; reflection);
                experiment_parameters(:interval, :QTCF; discrete=true)...)
            for (k, (method, result)) in enumerate(zip((:TCF, :QTCF), (tcf, qtcf)))
                record_discrete_run!(metadata, case.slug, method, result, grid,
                    grid.seed_node_counts[k], length(case.points))
            end
        end
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
    f, domain, algorithm = method === :QTCF ?
        (symmetric(conj; reflection=case.function_reflection),
            symmetric(domain; reflection=case.reflection), QTCF()) :
        (conj, domain, method === :AAA ? AAA() : TCF())
    approximate(f, domain, algorithm;
        experiment_parameters(:schwarz, method; discrete=points !== nothing)...)
end

function run_method(case, method, metadata; grid=nothing, seed_nodes=0)
    start = time_ns()
    result = fit_schwarz(case, method; points=grid === nothing ? nothing : grid.points)
    elapsed = (time_ns() - start) / 1e9
    rational = get_function(result)
    pole, residue = residues(rational)
    rotation = case.rotation
    pole .*= rotation
    residue .*= abs2(rotation)
    audit = curve_points(case.domain, SCHWARZ_AUDIT_N; periodic=RFA.isclosed(case.domain))
    error = max_error(z -> conj(rotation) * result(z / rotation), conj.(audit), audit)
    rational_degrees = degrees(rational)
    stop = status(result)
    grid === nothing || record_discrete_run!(metadata, case.slug, method, result, grid,
        seed_nodes, length(audit))
    (; pole, residue, error, elapsed, rational_degrees, failed=false,
        selected_nodes=length(nodes(rational)), recorded=length(result.history),
        returned_record=stop.best, sampled_error=stop.error, stop_reason=stop.reason)
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
    rows = Dict{Symbol,Any}()
    for (j, method) in enumerate(METHODS)
        print("  ", EXPERIMENT_MODE, " ", method, " ... ")
        flush(stdout)
        try
            row = run_method(case, method, metadata; grid,
                seed_nodes=grid === nothing ? 0 : grid.seed_node_counts[j])
            rows[method] = row
            @printf("error %.3e, degree %s, %.3f s\n",
                row.error, string(row.rational_degrees), row.elapsed)
        catch exception
            rows[method] = (
                pole=ComplexF64[], residue=ComplexF64[],
                error=NaN, elapsed=NaN, rational_degrees=(-1, -1), failed=true,
                selected_nodes=0, recorded=0, returned_record=0,
                sampled_error=NaN, stop_reason=:failed)
            println("failed: ", exception)
        end
        flush(stdout)
    end
    (; case, curve=curve_points(case.domain, SCHWARZ_CURVE_N), rows)
end

function run_experiments()
    global_logger(ConsoleLogger(stderr, Logging.Error))
    BLAS.set_num_threads(1)
    println("RationalFunctionApproximation ", pkgversion(RFA))
    println("Package source: ", pathof(RFA))
    println("Running all paper experiments.")
    EXPERIMENT_MODE === :discrete && println(
        "Discrete mode: one fixed shared grid per case, seeded by the continuum fits.")
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
