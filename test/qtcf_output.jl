# Figure, CSV, and HTML generation for both experiment modes.
using SHA, Dates, TOML

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
            cells = map(row) do value
                cell = string(value)
                occursin(r"[,\"\n]", cell) ? "\"" * replace(cell, "\"" => "\"\"") * "\"" : cell
            end
            println(io, join(cells, ','))
        end
    end
end

function record_discrete_run!(metadata, slug, method, result, grid, seed_nodes, validation_count)
    support = nodes(get_function(result))
    off_grid = count(z -> z ∉ grid.exact && !any(==(z), grid.points), support)
    off_grid == 0 || error("$slug/$method selected $off_grid off-grid nodes")
    stop = status(result)
    isequal(result.domain, grid.points) || error("$slug/$method changed the shared fitting grid")
    push!(metadata, (; case=slug, method,
        grid_source=grid.source,
        uniform_parameters=DISCRETE_UNIFORM_N,
        uniform_left_endpoint_included=true, uniform_right_endpoint_included=false,
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
        println(io, "domain,method,iteration,denominator_degree,max_error,best_validation,returned_record,sampled_error,stop_reason,returned_iteration,returned_sampled_error")
        for (column, panel) in enumerate(interval_results)
            axis = Axis(figure[1, column]; yscale=log10,
                title=LaTeXString(panel.title), xlabel="denominator degree",
                ylabel="max error", titlesize=14)
            for (run, color, label) in zip(panel.runs, colors, labels)
                data = history_data(run, panel.f, panel.points)
                stop = status(run)
                errors = max.(data.errors, floatmin(T))
                scatterlines!(axis, data.degrees, errors;
                    color, label, linewidth=1.0, markersize=2.8, alpha=0.65)
                scatter!(axis, [data.degrees[data.best]], [errors[data.best]];
                    color=:transparent, strokecolor=:black,
                    strokewidth=1.8, markersize=10)
                for j in eachindex(errors)
                    println(io, join((panel.slug, label, j,
                        data.degrees[j], data.errors[j], j == data.best,
                        j == stop.best, run.history[j].error, stop.reason,
                        stop.best, stop.error), ','))
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
        keep = isfinite.(row.pole) .& isfinite.(row.residue) .& .!iszero.(row.residue)
        append!(logs, log10.(abs.(row.residue[keep])))
    end
    isempty(logs) && return (-4.0, 0.0)
    lo = floor(Int, minimum(logs))
    hi = ceil(Int, maximum(logs))
    lo == hi && (lo -= 1)
    (Float64(lo), Float64(hi))
end

function residue_power_ticks(lo, hi)
    lower, upper = floor(Int, lo), ceil(Int, hi)
    step = max(1, cld(upper - lower, 6))
    exponents = collect(lower:step:upper)
    lo <= 0 <= hi && push!(exponents, 0)
    upper > 0 && push!(exponents, upper)
    unique!(sort!(exponents))
    Float64.(exponents), [LaTeXString("10^{$exponent}") for exponent in exponents]
end

function run_metadata()
    root = dirname(@__DIR__)
    git(args...) = try
        readchomp(Cmd(["git", "-C", root, args...]))
    catch
        "unavailable"
    end
    rows = Pair{String,String}[
        "mode" => string(EXPERIMENT_MODE),
        "generated_utc" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
        "julia_version" => string(VERSION),
        "float_type" => string(T), "float_precision_bits" => string(precision(T)),
        "blas_threads" => string(BLAS.get_num_threads()),
        "julia_threads" => string(Threads.nthreads()),
        "eval_method" => string(RFA.default_eval_method()),
        "weight_method" => string(RFA.default_weight_method()),
        "git_commit" => git("rev-parse", "HEAD"),
        "git_dirty" => string(!isempty(git("status", "--porcelain"))),
        "git_status" => git("status", "--porcelain"),
    ]
    manifest = joinpath(dirname(Base.active_project()), "Manifest.toml")
    if isfile(manifest)
        dependencies = TOML.parsefile(manifest)["deps"]
        for name in sort!(collect(keys(dependencies))), info in dependencies[name]
            push!(rows, "dependency/$name/$(info["uuid"])" =>
                "version=$(get(info, "version", ""));git-tree-sha1=$(get(info, "git-tree-sha1", ""))")
        end
    end
    for folder in ("src", "test"), (directory, _, files) in walkdir(joinpath(root, folder))
        for file in sort(files)
            (endswith(file, ".jl") || file in ("Project.toml", "Manifest.toml")) || continue
            path = joinpath(directory, file)
            push!(rows, "sha256/$(relpath(path, root))" => bytes2hex(sha256(read(path))))
        end
    end
    push!(rows, "sha256/Project.toml" => bytes2hex(sha256(read(joinpath(root, "Project.toml")))))
    write_csv("run_metadata.csv", "key,value", ((key, value) for (key, value) in rows))
end

function continuum_pole_figure(results, number;
    name="Figure_6_$(number)_continuum_poles_$(number - 1)", show_formula=true,
    colorrange, overview=false)
    lo, hi = colorrange

    nrows = length(results)
    row_step = show_formula ? 2 : 1
    figure = Figure(size=(720, 760), fontsize=8.5, figure_padding=2)
    for (j, method) in enumerate(METHODS)
        Label(figure[0, j], String(method); fontsize=13, font=:bold)
    end
    has_special = false
    for (i, result) in enumerate(results)
        case = result.case
        plot_row = row_step * (i - 1) + 1
        finite = [z for method in METHODS for z in result.rows[method].pole if isfinite(z)]
        extent = max(CONTINUUM_LIMITS[2], maximum(abs, real.(finite); init=0.0),
            maximum(abs, imag.(finite); init=0.0))
        limits = overview ? (-1.05extent, 1.05extent, -1.05extent, 1.05extent) : CONTINUUM_LIMITS
        for (j, method) in enumerate(METHODS)
            row = result.rows[method]
            axis = Axis(figure[plot_row, j]; aspect=DataAspect(), limits,
                title=CASE_LABELS[case.slug], titlesize=9,
                xlabel="Re z", ylabel="Im z", xlabelsize=8, ylabelsize=8,
                xticks=overview ? Makie.automatic : CONTINUUM_TICKS,
                yticks=overview ? Makie.automatic : CONTINUUM_TICKS,
                xticklabelsize=7, yticklabelsize=7)
            lines!(axis, real.(result.curve), imag.(result.curve);
                color=AAA_COLOR, linewidth=1.15)
            keep = isfinite.(row.pole) .& isfinite.(row.residue) .& .!iszero.(row.residue)
            if any(keep)
                scatter!(axis,
                    real.(row.pole[keep]), imag.(row.pole[keep]);
                    color=log10.(abs.(row.residue[keep])),
                    colormap=Reverse(:thermal), colorrange=(lo, hi), markersize=5.0)
            end
            zero_residue = isfinite.(row.pole) .& iszero.(row.residue)
            unknown = isfinite.(row.pole) .& .!isfinite.(row.residue)
            any(zero_residue) && scatter!(axis,
                real.(row.pole[zero_residue]), imag.(row.pole[zero_residue]);
                color=:transparent, strokecolor=:black, strokewidth=0.8, markersize=5.0)
            any(unknown) && scatter!(axis,
                real.(row.pole[unknown]), imag.(row.pole[unknown]);
                color=:gray40, marker=:xcross, markersize=6.0)
            outside = count(isfinite, row.pole) - count(in_window(row.pole, CONTINUUM_LIMITS))
            unavailable = count(unknown)
            unlocated = count(!isfinite, row.pole)
            has_special |= any(zero_residue) || unavailable > 0 || outside > 0 || unlocated > 0
            notes = String[]
            !overview && outside > 0 && push!(notes, "out: $outside")
            unavailable > 0 && push!(notes, "unavailable: $unavailable")
            unlocated > 0 && push!(notes, "nonfinite poles: $unlocated")
            !isempty(notes) && text!(axis, 0.02, 0.02; text=join(notes, "; "),
                space=:relative, align=(:left, :bottom), fontsize=6)
        end
        show_formula && Label(figure[2i, 1:3], LaTeXString(case.formula); fontsize=9)
    end
    tick_positions, tick_labels = residue_power_ticks(lo, hi)
    Colorbar(figure[1:row_step*nrows, 5]; limits=(lo, hi), colormap=Reverse(:thermal),
        label="residue", ticks=(tick_positions, tick_labels),
        labelsize=9, ticklabelsize=7)
    has_special && Label(figure[row_step*nrows+1, 1:3],
        "○ zero residue; × unavailable residue; out: poles shown in the overview"; fontsize=7)
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
    run_metadata()
    interval_outputs(interval_results)
    write_csv("schwarz_summary.csv",
        "case,method,max_error,degree,numerator_degree,denominator_degree,fit_seconds,failed,selected_nodes,recorded,returned_record,sampled_error,stop_reason,extraction_seconds,extraction_status,finite_poles,nonfinite_poles,nonfinite_residues",
        ((result.case.slug, method, row.error, row.rational_degrees[2],
            row.rational_degrees..., row.elapsed, row.failed, row.selected_nodes,
            row.recorded, row.returned_record, row.sampled_error, row.stop_reason,
            row.extraction_seconds, row.extraction_status, row.finite_poles,
            row.nonfinite_poles, row.nonfinite_residues)
            for result in schwarz_results for method in METHODS for row in (result.rows[method],)))
    cases = Tuple[(case.slug, :interval, (:TCF, :QTCF), length(case.points), 1.0 + 0.0im)
        for case in interval_cases]
    append!(cases, [(case.slug, :schwarz, METHODS, SCHWARZ_AUDIT_N, case.rotation) for case in qtcf_cases])
    stages = EXPERIMENT_MODE === :discrete ? (:continuum_seed, :discrete_fit) : (:continuum_fit,)
    write_csv("experiment_parameters.csv",
        "case,family,stage,method,float_type,tol,max_degree,max_iter,qtcf_initial_refinement,refinement,stagnation,allowed,loewner_tol,validation_points,blas_threads,coordinate_rotation,qtcf_start_parameter,qtcf_pairs_only",
        ((slug, family, stage, method, options.float_type, options.tol,
            get(options, :max_degree, ""), options.max_iter,
            get(options, :initial_refinement, ""), get(options, :refinement, ""),
            options.stagnation, options.allowed, get(options, :loewner_tol, ""),
            n, BLAS.get_num_threads(), rotation,
            method === :QTCF ? QTCF_START_PARAMETERS[slug] : "", method === :QTCF ? true : "")
            for (slug, family, methods, n, rotation) in cases for stage in stages for method in methods
            for options in (experiment_parameters(family, method; discrete=stage === :discrete_fit),)))
    if !isempty(metadata)
        write_csv("discrete_grid_metadata.csv", join(keys(first(metadata)), ','), values.(metadata))
    end
    for result in schwarz_results
        write_csv(result.case.slug * "_poles_residues.csv",
            "method,index,pole_real,pole_imag,residue_real,residue_imag,log10_abs_residue,in_plot_window,finite_pole,finite_residue",
            ((method, k, real(row.pole[k]), imag(row.pole[k]),
                real(row.residue[k]), imag(row.residue[k]),
                log10(abs(row.residue[k])), inside[k],
                isfinite(row.pole[k]), isfinite(row.residue[k]))
                for method in METHODS for row in (result.rows[method],)
                for inside in (in_window(row.pole, CONTINUUM_LIMITS),)
                for k in 1:min(length(row.pole), length(row.residue))))
    end
    colorrange = residue_colorrange(schwarz_results)
    overviews = Pair{String,String}[]
    for (i, first) in enumerate(1:4:length(schwarz_results))
        subset = schwarz_results[first:min(first + 3, end)]
        name = continuum_pole_figure(subset, i + 1; colorrange)
        if any(result -> any(method ->
            any(isfinite.(result.rows[method].pole) .&
                .!in_window(result.rows[method].pole, CONTINUUM_LIMITS)), METHODS), subset)
            overview = continuum_pole_figure(subset, i + 1;
                name=name * "_all_poles", colorrange, overview=true)
            push!(overviews, overview => "All finite poles: Schwarz curves $(first)–$(min(first+3, length(schwarz_results)))")
        else
            rm.(joinpath.(FIGURES, name * "_all_poles" .* (".pdf", ".png")); force=true)
        end
    end
    continuum_pole_figure(schwarz_results[2:2], 2;
        name="Figure_1", show_formula=false, colorrange)

    grid_link = isempty(metadata) ? "" :
        "<p><a href='numbers/discrete_grid_metadata.csv'>Shared fixed-grid metadata CSV</a></p>"
    figures = join(("<section><h2>$title</h2><p><a href='figures/$name.pdf'>PDF</a><a href='figures/$name.png'>PNG</a></p><img src='figures/$name.png' alt='$title'></section>"
        for (name, title) in (FIGURE_TITLES..., overviews...)), '\n')
    presets = join(("<li>$family: relative tolerance $(p.tol), degree limit $(p.max_degree), QTCF initial refinement $(p.initial_refinement), steady refinement $(p.refinement), stagnation $(p.stagnation).</li>"
        for (family, p) in pairs(EXPERIMENT_PARAMETERS)), '\n')
    write(joinpath(OUTPUT, "index.html"), """
        <!doctype html><html lang='en'><meta charset='utf-8'>
        <meta name='viewport' content='width=device-width,initial-scale=1'>
        <title>QTCF $EXPERIMENT_MODE experiments</title>
        <style>body{font:16px system-ui;max-width:1100px;margin:2em auto;padding:0 1em}
        img{max-width:100%;height:auto}section{margin:2em 0}a{margin-right:1em}</style>
        <h1>QTCF $EXPERIMENT_MODE experiments</h1>
        <p>All numerical examples; errors use the independent validation grids.</p>
        <p><a href='numbers/schwarz_summary.csv'>Schwarz summary CSV</a><a href='numbers/interval_convergence.csv'>Interval convergence CSV</a></p>
        <p>Each experiment family shares one $T preset across both modes, with unrestricted poles. Discrete fitting uses a fixed grid; iteration limits follow each method's degree convention.</p>
        <ul>$presets</ul>
        <p><a href='numbers/experiment_parameters.csv'>Applied fitting and grid-seed parameters</a><a href='numbers/run_metadata.csv'>Source and runtime metadata</a></p>
        <p>Residue colors use the full finite magnitude range shared by all panels. Hollow dots denote zero residues; gray crosses denote unavailable residues. Close-up panels show the stated window; the companion overviews include every finite pole with common limits across methods.</p>
        $grid_link
        $figures
        </html>
        """)
end
