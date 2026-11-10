# Assemble the Group L and P comparisons into a single PDF.
const CONTINUUM = @__DIR__
include("figure_base.jl")

const RELATIVE_TARGET = 4e-11
const L_GROUPS = (("ellipse_1_025", "ellipse_1_05", "upper_ellipse_1_05"),
    ("five_lobed", "superellipse_6", "analytic_inlet"))
const EXACT_GROUPS = (L_GROUPS..., INTERVALS, (S_CASES..., "odd_sine"))
const L_METHODS = (QTCF_MIDPOINT, PAIRED_CLASSIC)
const P_METHODS = (L_METHODS..., NATIVE_SQUARE)
const EXACT_METHODS = (L_METHODS, L_METHODS, P_METHODS, P_METHODS)
const GROUPED_FILENAME = "qtcf_required_comparison_figures.pdf"
const GROUP_P_FILENAME = "qtcf_required_comparison_group_p_figures.pdf"

function degree_histories(folder)
    degrees = Dict((row.case, row.method, row.block) => row
        for row in read_table(folder, "reduced_history_degrees"))
    histories = Dict{Tuple{String,String},Vector{NamedTuple}}()
    observations = read_table(folder, "histories")
    @assert length(observations) == length(degrees)
    for row in observations
        certified = degrees[(row.case, row.method, row.block)]
        @assert certified.status in ("certified_modular", "certified_exact")
        @assert row.q_degree == certified.nominal_degree && row.returned == certified.returned
        updated = merge(row, (; q_degree=certified.reduced_degree))
        push!(get!(Vector{NamedTuple}, histories, (row.case, row.method)), updated)
    end
    histories
end

function relative_histories(histories, summary)
    relative = Dict{Tuple{String,String},Vector{NamedTuple}}()
    for (key, rows) in histories
        scale = number(summary[key], :validation_target) / RELATIVE_TARGET
        @assert isfinite(scale) && scale > 0
        relative[key] = [merge(row, (;
            validation_error=string(number(row, :validation_error) / scale))) for row in rows]
    end
    relative
end

draw_target!(axis) = hlines!(axis, [RELATIVE_TARGET]; color=(:gray45, 0.70),
    linestyle=:dot, linewidth=0.9)

methods_with_history(histories, cases, methods) =
    [method for method in methods if any(case -> haskey(histories, (case, method)), cases)]

function pair_startup_note!(axis, data, histories, case, methods)
    if PAIRED_CLASSIC in methods && !haskey(histories, (case, PAIRED_CLASSIC))
        result = data.summary[(case, PAIRED_CLASSIC)]
        @assert result.returned_block == "0" && result.stop_reason == "nonfinite_weight"
        @assert result.failed_block == "1" && result.failed_leg == "2"
        text!(axis, 0.035, 0.085; text="$(STYLES[PAIRED_CLASSIC][1])\nEndpoint pair failed",
            space=:relative, align=(:left, :bottom), fontsize=8)
    end
end

function merge_pages(destination, pages)
    script = """
        import sys
        from pypdf import PdfWriter
        writer = PdfWriter()
        for path in sys.argv[2:]:
            writer.append(path, import_outline=False)
        writer.write(sys.argv[1])
        writer.close()
        """
    run(`python3 -c $script $destination $pages`)
end

function grouped_poles(data, exact_degrees, cases, destination; point_reflection=false)
    figure = Figure(size=(720, 760), fontsize=8.5, figure_padding=2)
    methods = point_reflection ? P_METHODS : L_METHODS
    columns = length(methods)
    panel_size = point_reflection ? 170 : 210
    for (j, method) in enumerate(methods)
        Label(figure[0, j], method_label(method);
            fontsize=point_reflection ? 11 : 12, font=:bold)
        colsize!(figure.layout, j, Fixed(panel_size))
    end
    has_special = false
    for (i, case) in enumerate(cases)
        plot_row = 2i - 1
        for (j, method) in enumerate(methods)
            axis = Axis(figure[plot_row, j]; aspect=DataAspect(),
                limits=(-POLE_WINDOW, POLE_WINDOW, -POLE_WINDOW, POLE_WINDOW),
                title=CASE_TITLE[case], titlesize=9, xlabel="Re z", ylabel="Im z",
                xlabelsize=8, ylabelsize=8, xticks=-2:2, yticks=-2:2,
                xticklabelsize=7, yticklabelsize=7)
            curve = data.domains[case]
            lines!(axis, real.(curve), imag.(curve); color=BLUE, linewidth=1.15)
            key = (case, method)
            candidates = get(data.poles, key, NamedTuple[])
            result = data.summary[key]
            exact_degree = get(exact_degrees, key, nothing)
            has_special |= pole_panel!(axis, candidates, result; exact_degree)
        end
        Label(figure[2i, 1:columns], FORMULA[case]; fontsize=9)
        rowsize!(figure.layout, plot_row, Fixed(panel_size))
        rowsize!(figure.layout, 2i, Fixed(17))
    end
    bottom = 2length(cases)
    Colorbar(figure[1:bottom, columns+2]; limits=(-4, 0), colormap=Reverse(:thermal), label="residue",
        ticks=(-4:0, [LaTeXString("10^{$k}") for k in -4:0]),
        labelsize=9, ticklabelsize=7, height=Relative(0.995), valign=:top)
    has_special && Label(figure[bottom+1, 1:columns], "○ zero residue    × unavailable residue"; fontsize=7)
    colsize!(figure.layout, columns+1, Fixed(24))
    rowgap!(figure.layout, 2)
    for i in 1:length(cases)-1
        rowgap!(figure.layout, 2i, 10)
    end
    colgap!(figure.layout, point_reflection ? 10 : 16)
    colgap!(figure.layout, columns, 0)
    colgap!(figure.layout, columns+1, 0)
    resize_to_layout!(figure)
    save(destination, figure)
end

function grouped_intervals(data, destination)
    figure = Figure(size=(720, 340), fontsize=11)
    highest_degree = maximum(parse(Int, row.q_degree)
        for case in INTERVALS for method in methods_with_history(data.histories, (case,), P_METHODS)
        for row in through_best_validation(data.histories[(case, method)]))
    xmax = max(100, 20 * cld(highest_degree, 20))
    flags = false
    for (column, case) in enumerate(INTERVALS)
        axis = Axis(figure[1, column]; yscale=log10, title=CASE_TITLE[case],
            xlabel="nominal denominator degree", ylabel="relative max error", titlesize=14, xticks=0:20:xmax)
        xlims!(axis, 0, xmax)
        draw_target!(axis)
        methods = methods_with_history(data.histories, (case,), P_METHODS)
        result = convergence_panel!(axis, data.histories, case, methods, :q_degree)
        flags |= result.nonfinite
        pair_startup_note!(axis, data, data.histories, case, P_METHODS)
        text!(axis, 0.96, 0.94; text=column == 1 ? "(a)" : "(b)", space=:relative,
            align=(:right, :top), fontsize=14)
    end
    methods = methods_with_history(data.histories, INTERVALS, P_METHODS)
    convergence_legend!(figure[2, 1:2], methods; crosses=flags)
    Label(figure[3, 1:2], "Circles mark the lowest independent-grid error. Dotted line marks relative tolerance 4 × 10⁻¹¹."; fontsize=9)
    rowgap!(figure.layout, 2)
    colgap!(figure.layout, 12)
    save(destination, figure)
end

function grouped_exact_page(data, histories, cases, methods, destination)
    figure = Figure(size=(850, 80 + 255 * length(cases)), fontsize=11)
    has_nonfinite = false
    for (i, case) in enumerate(cases)
        Label(figure[2i - 1, 1:2], CASE_TITLE[case]; fontsize=15, font=:bold)
        for (column, (field, label)) in enumerate(((:nodes, "Number of nodes"),
                (:q_degree, "exact reduced denominator degree")))
            axis = Axis(figure[2i, column]; yscale=log10, xlabel=label, ylabel="relative max error")
            draw_target!(axis)
            available = methods_with_history(histories, (case,), methods)
            flags = convergence_panel!(axis, histories, case, available, field)
            has_nonfinite |= flags.nonfinite
            column == 1 && pair_startup_note!(axis, data, histories, case, methods)
        end
    end
    legend_row = 2length(cases) + 1
    available = methods_with_history(histories, cases, methods)
    convergence_legend!(figure[legend_row, 1:2], available; crosses=has_nonfinite)
    note = "Circles mark the lowest independent-grid error. Dotted line marks relative tolerance 4 × 10⁻¹¹.\nCurves may overlap. Reduced degree can fall."
    has_nonfinite && (note *= "\nCrosses mark nonfinite errors on the independent grid.")
    Label(figure[legend_row + 1, 1:2], note; fontsize=9)
    rowgap!(figure.layout, 4)
    colgap!(figure.layout, 12)
    save(destination, figure)
end

function render_grouped_figures(; group_p_only=false)
    data = load_figures_data(CONTINUUM)
    histories = relative_histories(degree_histories(CONTINUUM), data.summary)
    exact_degrees = Dict(key => parse(Int, rows[returned_index(rows)].q_degree)
        for (key, rows) in histories)
    data = merge(data, (; histories=relative_histories(data.histories, data.summary)))
    mktempdir() do scratch
        pages = String[]
        if !group_p_only
            for (i, cases) in enumerate(L_GROUPS)
                path = joinpath(scratch, "line_$i.pdf")
                grouped_poles(data, exact_degrees, cases, path)
                push!(pages, path)
            end
        end
        point = joinpath(scratch, "P_poles.pdf")
        interval = joinpath(scratch, "P_intervals.pdf")
        grouped_poles(data, exact_degrees, S_CASES, point; point_reflection=true)
        grouped_intervals(data, interval)
        append!(pages, (point, interval))
        first_group = group_p_only ? length(L_GROUPS) + 1 : 1
        for i in first_group:length(EXACT_GROUPS)
            path = joinpath(scratch, "exact_$i.pdf")
            grouped_exact_page(data, histories, EXACT_GROUPS[i], EXACT_METHODS[i], path)
            push!(pages, path)
        end
        filename = group_p_only ? GROUP_P_FILENAME : GROUPED_FILENAME
        merge_pages(joinpath(@__DIR__, filename), pages)
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && render_grouped_figures(; group_p_only="--group-p-only" in ARGS)
