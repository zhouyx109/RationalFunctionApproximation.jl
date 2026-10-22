using RationalFunctionApproximation, ComplexRegions, CairoMakie, LaTeXStrings
using Logging, Printf, LinearAlgebra

isdefined(RFA, :QuadraticThiele) || error(
    "These experiments require the local QTCF extension. Run this file with --project=test from this checkout.")

const EXPERIMENT_MODE = isdefined(@__MODULE__, :QTCF_EXPERIMENT_MODE) ?
    QTCF_EXPERIMENT_MODE : :continuum
EXPERIMENT_MODE in (:continuum, :discrete) || error("Unknown experiment mode: $EXPERIMENT_MODE")

include("qtcf_parameters.jl")
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

    # Retain every distinct log sample and add each exact reflected coordinate.
    exact = Set(points)
    i = 1
    while i <= length(points)
        zr = coordinate(Z(reflection(points[i])))
        if zr ∉ exact
            push!(points, zr)
            push!(exact, zr)
        end
        i += 1
    end
    reflection_closed = all(z -> coordinate(Z(reflection(z))) in exact, points)
    reflection_closed || error("The shared fitting grid is not closed under reflection")
    println("  shared fixed grid: ", base_count, " base points, ", length(points),
        " after exact reflection closure (", length(unique(coordinate.(support))), " continuum support points)")
    (; points, exact, closed=RFA.isclosed(domain), source, base_count,
        reflection_closed, uniform_count=length(unique(coordinate.(uniform))),
        support_count=length(unique(coordinate.(support))), support_parameter_count=length(parameters),
        neighbor_count=length(unique(coordinate.(neighbors))),
        seed_node_counts=length.(nodes.(get_function.(seeds))))
end

function run_intervals(metadata)
    interval_seeds = map(interval_cases) do case
        println("\nInterval continuum seed: ", case.name)
        tcf = approximate(case.f, case.domain, TCF(); INTERVAL_TCF_PARAMETERS...)
        qtcf = approximate(symmetric(case.f), symmetric(case.domain);
            INTERVAL_QTCF_PARAMETERS...)
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
                discrete_parameters(INTERVAL_TCF_PARAMETERS)...)
            qtcf = approximate(symmetric(case.f), symmetric(grid.points; reflection);
                discrete_parameters(INTERVAL_QTCF_PARAMETERS)...)
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
    fit_domain = points === nothing ? case.domain : points
    rotation = EXPERIMENT_MODE === :continuum ? case.rotation : one(case.rotation)
    if rotation != 1
        # Use the same canonical coordinates for all three methods.
        original = fit_domain
        constructor = RFA.isclosed(original) ? ClosedCurve : Curve
        fit_domain = constructor(t -> point(original, t) / rotation,
            t -> ComplexRegions.tangent(original, t) / rotation, 0, 1)
    end
    # For Schwarz data, the value map must satisfy conj(ρ(z)) = τ(conj(z)).
    f, domain, algorithm = if method === :QTCF
        (symmetric(conj; reflection=rotation != 1 ? conj : case.function_reflection),
            symmetric(fit_domain; reflection=rotation != 1 ? conj : case.reflection), QTCF())
    else
        conj, fit_domain, method === :AAA ? AAA() : TCF()
    end
    parameters = SCHWARZ_PARAMETERS[method]
    points === nothing || (parameters = discrete_parameters(parameters))
    approximate(f, domain, algorithm; parameters...)
end

function run_method(case, method, metadata; grid=nothing, seed_nodes=0)
    start = time_ns()
    result = fit_schwarz(case, method; points=grid === nothing ? nothing : grid.points)
    elapsed = (time_ns() - start) / 1e9
    rational = get_function(result)
    pole, residue = residues(rational)
    rotation = EXPERIMENT_MODE === :continuum ? case.rotation : one(case.rotation)
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
        shared_discrete_grid(case.domain, case.reflection, seeds)
    end
    rows = Dict{Symbol,Any}()
    for (j, method) in enumerate(METHODS)
        print("  ", EXPERIMENT_MODE, " ", method, " ... ")
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
    println("Saved 4 figures as PDF and PNG, plus 11 CSV data files, to ", OUTPUT)
    (; interval_results, schwarz_results, metadata)
end

abspath(PROGRAM_FILE) == (@__FILE__) && run_experiments()
