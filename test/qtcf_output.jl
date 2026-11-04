# Figure, CSV, and HTML generation for both experiment modes.
using SHA, Dates, TOML

const OUTPUT = joinpath(@__DIR__,
    EXPERIMENT_MODE === :discrete ? "qtcf_discrete_output" : "qtcf_output")
const FIGURES = joinpath(OUTPUT, "figures")
const NUMBERS = joinpath(OUTPUT, "numbers")
const RUN_PROVENANCE = Ref{Any}(nothing)

utc_timestamp() = Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z"

function run_source_hashes(root, project)
    paths = Set{String}()
    for folder in ("src", "ext", "test")
        for (directory, directories, files) in walkdir(joinpath(root, folder))
            filter!(name -> name ∉ ("qtcf_output", "qtcf_discrete_output"), directories)
            for file in files
                (endswith(file, ".jl") || file == "Project.toml" ||
                    occursin(r"^Manifest(?:-v[0-9.]+)?\.toml$", file)) || continue
                push!(paths, joinpath(directory, file))
            end
        end
    end
    for directory in unique((root, dirname(project)))
        for file in readdir(directory)
            (file == "Project.toml" || occursin(r"^Manifest(?:-v[0-9.]+)?\.toml$", file)) || continue
            push!(paths, joinpath(directory, file))
        end
    end
    Dict("sha256/$(relpath(path, root))" => bytes2hex(sha256(read(path)))
        for path in sort!(collect(paths)))
end

function capture_run_provenance!()
    root = dirname(@__DIR__)
    project = Base.active_project()
    isnothing(project) && error("An active Julia project is required for experiment provenance.")
    git(args...) = try
        readchomp(Cmd(["git", "-C", root, args...]))
    catch
        "unavailable"
    end
    git_status = git("status", "--porcelain")
    rows = Pair{String,String}[
        "git_commit" => git("rev-parse", "HEAD"),
        "git_dirty" => git_status == "unavailable" ? "unavailable" : string(!isempty(git_status)),
        "git_status" => git_status,
        "active_project" => relpath(project, root),
    ]
    versioned = joinpath(dirname(project), "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
    manifest = isfile(versioned) ? versioned : joinpath(dirname(project), "Manifest.toml")
    if isfile(manifest)
        dependencies = get(TOML.parsefile(manifest), "deps", Dict())
        for name in sort!(collect(keys(dependencies))), info in dependencies[name]
            push!(rows, "dependency/$name/$(info["uuid"])" =>
                "version=$(get(info, "version", ""));git-tree-sha1=$(get(info, "git-tree-sha1", ""))")
        end
    end
    RUN_PROVENANCE[] = (; root, project, started_utc=utc_timestamp(),
        hashes=run_source_hashes(root, project), rows, source_unchanged=false,
        finished_utc="")
    nothing
end

function verify_run_provenance!()
    provenance = RUN_PROVENANCE[]
    isnothing(provenance) && error("Call capture_run_provenance!() before fitting experiments.")
    Base.active_project() == provenance.project || error("The active Julia project changed during the experiment run.")
    hashes = run_source_hashes(provenance.root, provenance.project)
    changed = sort!([path for path in union(keys(hashes), keys(provenance.hashes))
        if get(hashes, path, nothing) != get(provenance.hashes, path, nothing)])
    isempty(changed) || error("Experiment source or dependency files changed during fitting; outputs were not written: " * join(changed, ", "))
    RUN_PROVENANCE[] = merge(provenance, (; source_unchanged=true, finished_utc=utc_timestamp()))
    nothing
end

approximant_data(r) = r isa QuadraticThiele ? getfield(r, :data) : r

function save_figure(name, figure)
    save(joinpath(FIGURES, name * ".pdf"), figure)
    save(joinpath(FIGURES, name * ".png"), figure; px_per_unit=2)
end

function max_error(approximant, target, points)
    error = 0.0
    for i in eachindex(target, points)
        residual = abs(target[i] - approximant(points[i]))
        isfinite(residual) || return Inf
        error = max(error, residual)
    end
    return error
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
    off_grid = count(z -> sample_coordinate(z) ∉ grid.exact, support)
    off_grid == 0 || error("$slug/$method selected $off_grid off-grid nodes")
    stop = status(result)
    isequal(result.domain, grid.points) || error("$slug/$method changed the shared fitting grid")
    push!(metadata, (; case=slug, method,
        grid_source=grid.source,
        uniform_parameters=iszero(grid.uniform_count) ? 0 : get(grid, :uniform_subdivisions, DISCRETE_UNIFORM_N),
        uniform_left_endpoint_included=!iszero(grid.uniform_count),
        uniform_right_endpoint_included=!iszero(grid.uniform_count) && !grid.closed,
        uniform_points=grid.uniform_count,
        continuum_support_points=grid.support_count, seed_nodes,
        continuum_node_scope=iszero(seed_nodes) ? "none" : "all_recorded_stages",
        support_parameters=grid.support_parameter_count,
        log_exponents=iszero(grid.neighbor_count) ? "" : join(DISCRETE_LOG_EXPONENTS, ";"),
        neighbor_scales=iszero(grid.neighbor_count) ? 0 : length(DISCRETE_LOG_EXPONENTS),
        newman_n=get(grid, :newman_n, 0), newman_base_n=get(grid, :newman_base_n, 0),
        newman_density=get(grid, :newman_density, 0),
        newman_base_eta=get(grid, :newman_base_eta, NaN), newman_eta=get(grid, :newman_eta, NaN),
        newman_minimum_positive=get(grid, :newman_minimum_positive, NaN),
        neighbor_points=grid.neighbor_count, base_points=grid.base_count,
        shared_points=length(grid.points), reflection_added=length(grid.points)-grid.base_count,
        closed_curve=grid.closed, reflection_closed=grid.reflection_closed,
        initial_node_real=real(first(support)), initial_node_imag=imag(first(support)),
        selected_nodes=length(support), off_grid_nodes=off_grid, validation_points=validation_count,
        stop_reason=stop.reason, iterations=length(result.history), best_iteration=stop.best))
end

# Fixed validation grids shared by all methods, records, and fitting modes.
function validation_grid(case)
    points = hasproperty(case, :points) ? case.points :
        curve_points(case.fit_domain, SCHWARZ_AUDIT_N; periodic=RFA.isclosed(case.fit_domain))
    checksum = bytes2hex(sha256(reinterpret(UInt8, vec(ComplexF64.(points)))))
    (; points, metadata=(grid_sha256=checksum,))
end

function history_data(result, points, target)
    denominator_degrees = [last(degrees(record.interpolant)) for record in result.history]
    errors = [max_error(record.interpolant, target, points) for record in result.history]
    error, best = findmin(errors)
    isfinite(error) || error("No finite interval approximation was found.")
    (; degrees=denominator_degrees, errors, best)
end

function interval_outputs(interval_results)
    figure = Figure(size=(720, 340), fontsize=11)
    colors = (TCF_COLOR, QTCF_COLOR)
    labels = ("TCF", "symmetric QTCF")

    open(joinpath(NUMBERS, "interval_convergence.csv"), "w") do io
        println(io, "domain,method,iteration,denominator_degree,max_error,best_validation,returned_record,sampled_error,stop_reason,returned_iteration,returned_sampled_error,representation,plotted,validation_points,validation_nodes,validation_grid_sha256")
        for (column, panel) in enumerate(interval_results)
            grid = validation_grid(panel)
            target = panel.f.(grid.points)
            axis = Axis(figure[1, column]; yscale=log10,
                title=LaTeXString(panel.title), xlabel="denominator degree",
                ylabel="max error", titlesize=14,
                xticks=0:20:INTERVAL_DISPLAY_DEGREE)
            xlims!(axis, 0, INTERVAL_DISPLAY_DEGREE)
            for (run, color, label) in zip(panel.runs, colors, labels)
                data = history_data(run, grid.points, target)
                stop = status(run)
                prefix = 1:data.best
                errors = max.(data.errors[prefix], floatmin(T))
                finite_errors = ifelse.(isfinite.(errors), errors, NaN)
                scatterlines!(axis, data.degrees[prefix], finite_errors;
                    color, label, linewidth=1.0, markersize=2.8, alpha=0.65)
                nonfinite = findall(!isfinite, data.errors[prefix])
                if !isempty(nonfinite)
                    ceiling = 10maximum(filter(isfinite, errors); init=one(T))
                    scatter!(axis, data.degrees[nonfinite], fill(ceiling, length(nonfinite));
                        color, marker=:xcross, markersize=7)
                    text!(axis, 0.02, 0.98; text="× nonfinite error", space=:relative,
                        align=(:left, :top), fontsize=9)
                end
                scatter!(axis, [data.degrees[data.best]], [errors[end]];
                    color=:transparent, strokecolor=:black,
                    strokewidth=1.8, markersize=10)
                for j in eachindex(data.errors)
                    println(io, join((panel.slug, label, j,
                        data.degrees[j], data.errors[j], j == data.best,
                        j == stop.best, run.history[j].error, stop.reason,
                        stop.best, stop.error,
                        nameof(typeof(approximant_data(run.history[j].interpolant))),
                        j <= data.best && data.degrees[j] <= INTERVAL_DISPLAY_DEGREE,
                        length(grid.points), 0, grid.metadata.grid_sha256), ','))
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
    provenance = RUN_PROVENANCE[]
    (isnothing(provenance) || !provenance.source_unchanged) &&
        error("Experiment provenance must be verified before writing metadata.")
    rows = Pair{String,String}[
        "mode" => string(EXPERIMENT_MODE),
        "started_utc" => provenance.started_utc,
        "finished_utc" => provenance.finished_utc,
        "generated_utc" => utc_timestamp(),
        "source_unchanged" => string(provenance.source_unchanged),
        "julia_version" => string(VERSION),
        "float_type" => string(T), "float_precision_bits" => string(precision(T)),
        "blas_threads" => string(BLAS.get_num_threads()),
        "julia_threads" => string(Threads.nthreads()),
        "eval_method" => string(RFA.default_eval_method()),
        "weight_method/AAA" => "not applicable",
        "weight_method/TCF" => string(EXPERIMENT_WEIGHT_METHODS.TCF),
        "weight_method/QTCF" => string(EXPERIMENT_WEIGHT_METHODS.QTCF),
        "sample_evaluation/AAA" => "native barycentric evaluation",
        "sample_evaluation/TCF" => "native array evaluation",
        "sample_evaluation/QTCF" =>
            "unscaled array sweep with scalar NaN retries and selected-node common-zero recovery",
        "validation_evaluation" => "scalar broadcast on fixed paper test sets",
        "validation/interval" => "X1 union X2 union (-X2) union (X2-1); imaginary uses iX",
        "validation/interval_uniform_points" => string(ABS_UNIFORM_N),
        "validation/interval_log_exponents" => join(ABS_LOG_EXPONENTS, ";"),
        "validation/schwarz_uniform_points" => string(SCHWARZ_AUDIT_N),
    ]
    append!(rows, provenance.rows)
    append!(rows, sort!(collect(provenance.hashes); by=first))
    write_csv("run_metadata.csv", "key,value", ((key, value) for (key, value) in rows))
end

function continuum_pole_figure(results, name; show_formula=true)
    lo, hi = POLE_RESIDUE_LOG_RANGE

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
        for (j, method) in enumerate(METHODS)
            row = result.rows[method]
            axis = Axis(figure[plot_row, j]; aspect=DataAspect(), limits=CONTINUUM_LIMITS,
                title=CASE_LABELS[case.slug], titlesize=9,
                xlabel="Re z", ylabel="Im z", xlabelsize=8, ylabelsize=8,
                xticks=CONTINUUM_TICKS, yticks=CONTINUUM_TICKS,
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
            unavailable = count(unknown)
            unlocated = count(!isfinite, row.pole)
            has_special |= any(zero_residue) || unavailable > 0 || unlocated > 0
            notes = String[]
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
        labelsize=9, ticklabelsize=7,
        height=Relative(nrows > 1 ? 0.995 : 1.0), valign=:top)
    has_special && Label(figure[row_step*nrows+1, 1:3],
        "○ zero residue; × unavailable residue"; fontsize=7)
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

function pole_figures(schwarz_results)
    for (i, first) in enumerate(1:4:length(schwarz_results))
        subset = schwarz_results[first:min(first + 3, end)]
        continuum_pole_figure(subset, "Figure_6_$(i + 1)_continuum_poles_$i")
    end
    ellipse = only(result for result in schwarz_results if result.case.slug == "ellipse_1_05")
    continuum_pole_figure((ellipse,), "Figure_1"; show_formula=false)
end

function write_outputs(interval_results, schwarz_results, metadata)
    verify_run_provenance!()
    mkpath.((FIGURES, NUMBERS))
    run_metadata()
    interval_outputs(interval_results)
    write_csv("schwarz_summary.csv",
        "case,method,max_error,degree,numerator_degree,denominator_degree,fit_seconds,failed,selected_nodes,recorded,returned_record,sampled_error,stop_reason,extraction_seconds,extraction_status,finite_poles,nonfinite_poles,nonfinite_residues,validation_points,validation_nodes,validation_grid_sha256",
        ((result.case.slug, method, row.error, row.rational_degrees[2],
            row.rational_degrees..., row.elapsed, row.failed, row.selected_nodes,
            row.recorded, row.returned_record, row.sampled_error, row.stop_reason,
            row.extraction_seconds, row.extraction_status, row.finite_poles,
            row.nonfinite_poles, row.nonfinite_residues,
            row.validation_points, row.validation_nodes, row.validation_grid_sha256)
            for result in schwarz_results for method in METHODS for row in (result.rows[method],)))
    cases = Tuple[(case.slug, :interval, (:TCF, :QTCF), 1.0 + 0.0im)
        for case in interval_cases]
    append!(cases, [(case.slug, :schwarz, METHODS, case.rotation) for case in qtcf_cases])
    stages = EXPERIMENT_MODE === :discrete ? (:continuum_seed, :discrete_fit) : (:continuum_fit,)
    write_csv("experiment_parameters.csv",
        "case,family,stage,method,float_type,tol,degree_limit,max_degree,max_iter,qtcf_initial_refinement,refinement,stagnation,allowed,loewner_tol,validation_points,blas_threads,coordinate_rotation,qtcf_start_parameters,qtcf_pairs_only,validation_rule,validation_log_exponents",
        ((slug, family, stage, method, options.float_type, options.tol,
            EXPERIMENT_PARAMETERS[family].max_degree,
            get(options, :max_degree, ""), options.max_iter,
            get(options, :initial_refinement, ""), get(options, :refinement, ""),
            options.stagnation, options.allowed, get(options, :loewner_tol, ""),
            family === :interval ? length(abs_check) : SCHWARZ_AUDIT_N,
            BLAS.get_num_threads(), rotation,
            method === :QTCF ? repr(startup) : "", method === :QTCF ? true : "",
            family === :interval ? "X1+X2+(-X2)+(X2-1)" : "uniform_curve_parameter",
            family === :interval ? join(ABS_LOG_EXPONENTS, ";") : "")
            for (slug, family, methods, rotation) in cases
            for stage in (EXPERIMENT_MODE === :discrete && family === :interval ? (:discrete_fit,) : stages)
            for method in methods
            for options in (experiment_parameters(family, method; discrete=stage === :discrete_fit),)
            for startup in ((QTCF_START_PARAMETER,),)))
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
                for k in eachindex(row.pole, row.residue)))
    end
    pole_figures(schwarz_results)

    grid_link = isempty(metadata) ? "" :
        "<p>Discrete Schwarz grids use $DISCRETE_UNIFORM_N uniform subdivisions, continuum interpolation nodes, and neighbors at parameter offsets 10<sup>k</sup> with k in $DISCRETE_LOG_EXPONENTS. <a href='numbers/discrete_grid_metadata.csv'>Shared fixed-grid metadata CSV</a></p>"
    figures = join(("<section><h2>$title</h2><p><a href='figures/$name.pdf'>PDF</a><a href='figures/$name.png'>PNG</a></p><img src='figures/$name.png' alt='$title'></section>"
        for (name, title) in FIGURE_TITLES), '\n')
    presets = join(("<li>$family, $method: relative tolerance $(options.tol), degree limit $(p.max_degree), initial refinement $(get(options, :initial_refinement, 15)), steady refinement $(options.refinement), stagnation $(options.stagnation), iteration limit $(options.max_iter).</li>"
        for (family, p) in pairs(EXPERIMENT_PARAMETERS)
        for method in (family === :interval ? (:TCF, :QTCF) : METHODS)
        for options in (experiment_parameters(family, method),)), '\n')
    write(joinpath(OUTPUT, "index.html"), """
        <!doctype html><html lang='en'><meta charset='utf-8'>
        <meta name='viewport' content='width=device-width,initial-scale=1'>
        <title>QTCF $EXPERIMENT_MODE experiments</title>
        <style>body{font:16px system-ui;max-width:1100px;margin:2em auto;padding:0 1em}
        img{max-width:100%;height:auto}section{margin:2em 0}a{margin-right:1em}</style>
        <h1>QTCF $EXPERIMENT_MODE experiments</h1>
        <p>Interval errors use the paper's fixed set X = X₁ ∪ X₂ ∪ (−X₂) ∪ (X₂−1), where X₁ has 10,001 uniform points and X₂ = {2<sup>−0.1k</sup>: k = 10,…,1000}. The imaginary interval uses iX.</p>
        <p>Schwarz errors use $SCHWARZ_AUDIT_N uniform parameter points on each curve. Closed curves count their common endpoint once. Validation points are independent of the approximant and shared by all methods and both fitting modes.</p>
        <p><a href='numbers/schwarz_summary.csv'>Schwarz summary CSV</a><a href='numbers/interval_convergence.csv'>Interval convergence CSV</a></p>
        <p>Each method uses the same $T tolerance and iteration budget in both sampling modes, with unrestricted poles. Refinement applies only to continuum sampling. Discrete intervals preserve the original TCF fitting grid and union order; Schwarz grids include continuum interpolation nodes.</p>
        <p>QTCF starts with one scalar block at the parameter midpoint, followed by symmetric pairs.</p>
        <ul>$presets</ul>
        <p><a href='numbers/experiment_parameters.csv'>Applied fitting and grid-seed parameters</a><a href='numbers/run_metadata.csv'>Source and runtime metadata</a></p>
        <p>Interval curves end at their minimum validation error, marked by an outlined point, and display denominator degrees 0–$(INTERVAL_DISPLAY_DEGREE). All history records remain in the CSVs. Pole panels use the paper's fixed window and residue colors from 10<sup>-4</sup> to 1, with colors saturated outside this range. Hollow dots denote zero residues; gray crosses denote unavailable residues. Every returned pole and residue remains in the CSVs.</p>
        $grid_link
        $figures
        </html>
        """)
end
