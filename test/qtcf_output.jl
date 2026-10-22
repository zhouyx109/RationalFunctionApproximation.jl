# Figure, CSV, and HTML generation for both experiment modes.
const OUTPUT = joinpath(@__DIR__,
    EXPERIMENT_MODE === :discrete ? "qtcf_discrete_output" : "qtcf_output")
const FIGURES = joinpath(OUTPUT, "figures")
const NUMBERS = joinpath(OUTPUT, "numbers")

function save_figure(name, figure)
    save(joinpath(FIGURES, name * ".pdf"), figure)
    save(joinpath(FIGURES, name * ".png"), figure; px_per_unit=2)
end

function max_error(approximant, target, points)
    errors = abs.(target .- approximant.(points))
    return all(isfinite, errors) ? maximum(errors; init=0.0) : Inf
end

function write_csv(name, header, rows)
    open(joinpath(NUMBERS, name), "w") do io
        println(io, header)
        for row in rows
            println(io, join(row, ','))
        end
    end
end

function record_discrete_run!(metadata, slug, method, result, grid, seed_nodes, validation_count)
    support = nodes(get_function(result))
    off_grid = count(z -> z ∉ grid.exact && !any(==(z), grid.points), support)
    off_grid == 0 || error("$slug/$method selected $off_grid off-grid nodes")
    stop = status(result)
    push!(metadata, (; case=slug, method,
        grid_source=grid.source,
        uniform_parameters=DISCRETE_UNIFORM_N, uniform_endpoint_included=false,
        uniform_points=grid.uniform_count,
        continuum_support_points=grid.support_count, seed_nodes,
        support_parameters=grid.support_parameter_count,
        log_min=DISCRETE_LOG_MIN, log_step=DISCRETE_LOG_STEP, log_max=DISCRETE_LOG_MAX,
        neighbor_points=grid.neighbor_count, base_points=grid.base_count,
        shared_points=length(grid.points), reflection_added=length(grid.points)-grid.base_count,
        closed_curve=grid.closed, reflection_closed=grid.reflection_closed,
        selected_nodes=length(support), off_grid_nodes=off_grid, validation_points=validation_count,
        stop_reason=stop.reason, iterations=length(result.history), best_iteration=stop.best))
end

function history_data(result, f, points)
    target = f.(points)
    approximants = getproperty.(result.history, :interpolant)
    denominator_degrees = [degrees(r)[2] for r in approximants]
    errors = Float64[max_error(r, target, points) for r in approximants]
    candidates = findall(isfinite, errors)
    isempty(candidates) && error("No finite interval approximation was found.")
    best = candidates[argmin(errors[candidates])]
    (; degrees=denominator_degrees, errors, best)
end

function interval_outputs(interval_results)
    figure = Figure(size=(720, 340), fontsize=11)
    colors = (TCF_COLOR, QTCF_COLOR)
    labels = ("TCF", "symmetric QTCF")

    open(joinpath(NUMBERS, "interval_convergence.csv"), "w") do io
        println(io, "domain,method,iteration,denominator_degree,max_error,best")
        for (column, panel) in enumerate(interval_results)
            axis = Axis(figure[1, column]; yscale=log10,
                title=LaTeXString(panel.title), xlabel="denominator degree",
                ylabel="max error", titlesize=14)
            for (run, color, label) in zip(panel.runs, colors, labels)
                data = history_data(run, panel.f, panel.points)
                errors = max.(data.errors, floatmin(T))
                scatterlines!(axis, data.degrees, errors;
                    color, label, linewidth=1.0, markersize=2.8, alpha=0.65)
                scatter!(axis, [data.degrees[data.best]], [errors[data.best]];
                    color=:transparent, strokecolor=:black,
                    strokewidth=1.8, markersize=10)
                for j in eachindex(errors)
                    println(io, join((panel.slug, label, j,
                        data.degrees[j], data.errors[j], j == data.best), ','))
                end
            end
            text!(axis, 0.96, 0.94; text=column == 1 ? "(a)" : "(b)",
                space=:relative, align=(:right, :top), fontsize=14)
        end
    end

    legend_elements = [
        [LineElement(color=color, linewidth=1.3),
            MarkerElement(color=color, marker=:circle, markersize=4)]
        for color in colors
    ]
    Legend(figure[2, 1:2], legend_elements, collect(labels);
        orientation=:horizontal, framevisible=true, padding=(4, 4, 2, 2),
        labelsize=10)
    rowgap!(figure.layout, 2)
    colgap!(figure.layout, 12)
    save_figure("Figure_6_1_interval_convergence", figure)
end

function in_window(poles, limits)
    xmin, xmax, ymin, ymax = limits
    isfinite.(poles) .&
        (xmin .<= real.(poles)) .& (real.(poles) .<= xmax) .&
        (ymin .<= imag.(poles)) .& (imag.(poles) .<= ymax)
end

function residue_colorrange(schwarz_results)
    logs = Float64[]
    for result in schwarz_results, method in METHODS
        row = result.rows[method]
        keep = in_window(row.pole, CONTINUUM_LIMITS) .& isfinite.(row.residue)
        append!(logs, log10.(abs.(row.residue[keep]) .+ eps(T)))
    end
    # Use shared residue limits so every panel has the same color scale.
    sort!(logs)
    if isempty(logs)
        lo, hi = -4, 0
    else
        n = length(logs)
        lo = floor(Int, logs[clamp(ceil(Int, 0.10n), 1, n)])
        hi = ceil(Int, logs[clamp(ceil(Int, 0.99n), 1, n)])
        lo == hi && (lo -= 1)
    end
    (Float64(lo), Float64(hi))
end

function residue_power_ticks(lo, hi)
    lower = floor(Int, lo)
    upper = ceil(Int, hi)
    step = max(1, cld(abs(lower), 4))
    exponents = collect(lower:step:-1)
    lo <= 0 <= hi && push!(exponents, 0)
    upper > 0 && push!(exponents, upper)
    unique!(sort!(exponents))
    Float64.(exponents), [LaTeXString("10^{$exponent}") for exponent in exponents]
end

function continuum_pole_figure(results, number;
    name="Figure_6_$(number)_continuum_poles_$(number - 1)", show_formula=true, colorrange)
    lo, hi = colorrange

    nrows = length(results)
    row_step = show_formula ? 2 : 1
    figure = Figure(size=(720, 760), fontsize=8.5, figure_padding=2)
    for (j, method) in enumerate(METHODS)
        Label(figure[0, j], String(method); fontsize=13, font=:bold)
    end
    last_scatter = nothing
    for (i, result) in enumerate(results)
        case = result.case
        plot_row = row_step * (i - 1) + 1
        for (j, method) in enumerate(METHODS)
            row = result.rows[method]
            axis = Axis(figure[plot_row, j]; aspect=DataAspect(), limits=CONTINUUM_LIMITS,
                title=CASE_LABELS[case.slug], titlesize=9,
                xlabel="Re z", ylabel="Im z", xlabelsize=8, ylabelsize=8,
                xticks=CONTINUUM_TICKS, yticks=CONTINUUM_TICKS,
                xticklabelsize=7, yticklabelsize=7)
            lines!(axis, real.(result.curve), imag.(result.curve);
                color=AAA_COLOR, linewidth=1.15)
            keep = in_window(row.pole, CONTINUUM_LIMITS) .& isfinite.(row.residue)
            if any(keep)
                last_scatter = scatter!(axis,
                    real.(row.pole[keep]), imag.(row.pole[keep]);
                    color=log10.(abs.(row.residue[keep]) .+ eps(T)),
                    colormap=Reverse(:thermal), colorrange=(lo, hi), markersize=5.0)
            end
        end
        show_formula && Label(figure[2i, 1:3], LaTeXString(case.formula); fontsize=9)
    end
    tick_positions, tick_labels = residue_power_ticks(lo, hi)
    last_scatter !== nothing && Colorbar(figure[1:row_step*nrows, 5], last_scatter;
        label="residue", ticks=(tick_positions, tick_labels),
        labelsize=9, ticklabelsize=7)
    for column in 1:3
        colsize!(figure.layout, column, Fixed(175))
    end
    colsize!(figure.layout, 4, Fixed(24))
    for row in 1:nrows
        rowsize!(figure.layout, row_step * (row - 1) + 1, Fixed(175))
        show_formula && rowsize!(figure.layout, 2row, Fixed(17))
    end
    rowgap!(figure.layout, 2)
    for row in 1:nrows-1
        rowgap!(figure.layout, row_step * row, 10)
    end
    colgap!(figure.layout, 16)
    colgap!(figure.layout, 3, 0)
    colgap!(figure.layout, 4, 0)
    resize_to_layout!(figure)

    save_figure(name, figure)
    name
end

function write_outputs(interval_results, schwarz_results, metadata)
    mkpath.((FIGURES, NUMBERS))
    interval_outputs(interval_results)
    write_csv("schwarz_summary.csv",
        "case,method,max_error,degree,numerator_degree,denominator_degree,time_seconds,failed,selected_nodes,recorded,returned_record,sampled_error,stop_reason",
        ((result.case.slug, method, row.error, row.rational_degrees[2],
            row.rational_degrees..., row.elapsed, row.failed, row.selected_nodes,
            row.recorded, row.returned_record, row.sampled_error, row.stop_reason)
            for result in schwarz_results for method in METHODS for row in (result.rows[method],)))
    if EXPERIMENT_MODE === :continuum
        cases = [(case.slug, [:TCF, :QTCF], length(case.points), 1.0 + 0.0im)
            for case in interval_cases]
        append!(cases, [(case.slug, collect(METHODS), SCHWARZ_AUDIT_N, case.rotation) for case in qtcf_cases])
        write_csv("experiment_parameters.csv",
            "case,method,float_type,tol,max_degree,max_iter,initial_refinement,refinement,stagnation,allowed,validation_points,blas_threads,coordinate_rotation",
            ((slug, method, T, p.tol, p.max_degree,
                method === :TCF ? 2p.max_degree : p.max_degree,
                method === :QTCF ? p.initial_refinement : RFA_INITIAL_REFINEMENT,
                p.refinement, p.stagnation, p.allowed, n, BLAS.get_num_threads(), rotation)
                for (slug, methods, n, rotation) in cases for method in methods
                for p in (method === :QTCF ? CONTINUUM_QTCF_PARAMETERS : CONTINUUM_PARAMETERS,)))
    elseif !isempty(metadata)
        write_csv("discrete_grid_metadata.csv", join(keys(first(metadata)), ','), values.(metadata))
    end
    for result in schwarz_results
        write_csv(result.case.slug * "_poles_residues.csv",
            "method,index,pole_real,pole_imag,residue_real,residue_imag,log10_abs_residue,in_plot_window",
            ((method, k, real(row.pole[k]), imag(row.pole[k]),
                real(row.residue[k]), imag(row.residue[k]),
                log10(abs(row.residue[k]) + eps(T)), inside[k])
                for method in METHODS for row in (result.rows[method],)
                for inside in (in_window(row.pole, CONTINUUM_LIMITS),)
                for k in 1:min(length(row.pole), length(row.residue))))
    end
    colorrange = residue_colorrange(schwarz_results)
    for (i, first) in enumerate(1:4:length(schwarz_results))
        continuum_pole_figure(schwarz_results[first:min(first + 3, end)], i + 1; colorrange)
    end
    continuum_pole_figure(schwarz_results[2:2], 2;
        name="Figure_1", show_formula=false, colorrange)

    open(joinpath(OUTPUT, "index.html"), "w") do io
        println(io, "<!doctype html><html lang='en'><meta charset='utf-8'>",
            "<meta name='viewport' content='width=device-width,initial-scale=1'>",
            "<title>QTCF ", EXPERIMENT_MODE, " experiments</title>",
            "<style>body{font:16px system-ui;max-width:1100px;margin:2em auto;padding:0 1em}",
            "img{max-width:100%;height:auto}section{margin:2em 0}a{margin-right:1em}</style>",
            "<h1>QTCF ", EXPERIMENT_MODE, " experiments</h1>",
            "<p>All numerical examples; errors use the independent validation grids.</p>",
            "<p><a href='numbers/schwarz_summary.csv'>Schwarz summary CSV</a>",
            "<a href='numbers/interval_convergence.csv'>Interval convergence CSV</a></p>")
        p = CONTINUUM_PARAMETERS
        EXPERIMENT_MODE === :continuum && println(io,
            "<p>Every continuum case uses ", T, ", relative tolerance ", p.tol,
            ", degree limit ", p.max_degree,
            ", initial refinement ", CONTINUUM_QTCF_PARAMETERS.initial_refinement,
            ", steady refinement ", p.refinement, ", stagnation ", p.stagnation, ", ",
            "and unrestricted poles. Iteration limits follow each method's degree convention.</p>",
            "<p><a href='numbers/experiment_parameters.csv'>Complete shared parameter preset</a></p>")
        EXPERIMENT_MODE === :discrete && println(io,
            "<p><a href='numbers/discrete_grid_metadata.csv'>Shared fixed-grid metadata CSV</a></p>")
        for (name, title) in FIGURE_TITLES
            println(io, "<section><h2>", title, "</h2><p><a href='figures/", name,
                ".pdf'>PDF</a><a href='figures/", name, ".png'>PNG</a></p><img src='figures/",
                name, ".png' alt='", title, "'></section>")
        end
        println(io, "</html>")
    end
end
