# Read the experiment tables and provide shared styles for the error and pole plots.
using CairoMakie, LaTeXStrings, DelimitedFiles

# Match the CairoMakie styling in 03_QTCF/test/qtcf_output.jl.
const INTERVALS = ("real", "imaginary")
const S_CASES = ("smooth_s", "two_semicircle_s")
const POLE_WINDOW = 2.25
const CASE_TITLE = Dict(
    "real" => L"|x|\ \mathrm{on}\ [-1,1]",
    "imaginary" => L"|z|\ \mathrm{on}\ [-i,i]",
    "ellipse_1_025" => "Ellipse (1, 1/4)",
    "ellipse_1_05" => "Ellipse (1, 1/2)",
    "upper_ellipse_1_05" => "Upper ellipse (1, 1/2)",
    "five_lobed" => "Five-lobed curve",
    "smooth_s" => "Smooth S",
    "two_semicircle_s" => "Two-semicircle S",
    "superellipse_6" => "Superellipse",
    "analytic_inlet" => "Analytic inlet",
    "odd_sine" => "Smooth odd real target",
)
const FORMULA = Dict(
    "ellipse_1_025" => L"z(t)=\cos t+\frac{i}{4}\sin t",
    "ellipse_1_05" => L"z(t)=\cos t+\frac{i}{2}\sin t",
    "upper_ellipse_1_05" => L"z(t)=\cos(\pi(1-t))+\frac{i}{2}\sin(\pi(1-t))",
    "five_lobed" => L"z(t)=(1+0.2\sin 5t)e^{it}",
    "smooth_s" => L"z(t)=-0.75\sin(2\pi t)+i(2-4t)",
    "two_semicircle_s" => L"z(t)=i+e^{i(\pi/2+2\pi t)}\ (t\leq1/2),\quad z(t)=-i+e^{i(3\pi/2-2\pi t)}\ (t>1/2)",
    "superellipse_6" => L"z(t)=e^{2\pi it}/(\cos^6(2\pi t)+\sin^6(2\pi t))^{1/6}",
    "analytic_inlet" => L"z(t)=0.1-0.1\cos t-0.7\cos2t-0.1\cos3t-i(0.8\sin t+0.4\sin2t+0.2\sin3t)",
)
const BLUE = RGBf(0 / 255, 114 / 255, 178 / 255)
const GREEN = RGBf(0 / 255, 158 / 255, 115 / 255)
const VERMILION = RGBf(213 / 255, 94 / 255, 0 / 255)
const SKY = RGBf(86 / 255, 180 / 255, 233 / 255)
const QTCF_MIDPOINT = "QTCF"
const PAIRED_CLASSIC = "forced_Classic"
const NATIVE_SQUARE = "native_squared"
const STYLES = Dict(
    QTCF_MIDPOINT => ("QTCF (midpoint)", GREEN, :solid, :circle),
    PAIRED_CLASSIC => ("Paired TCF", VERMILION, :solid, :utriangle),
    NATIVE_SQUARE => ("Native square TCF", SKY, :dashdot, :diamond),
)

function read_table(folder, name)
    rows, header = readdlm(joinpath(folder, name * ".csv"), ',', String, '\n'; header=true)
    keys = Tuple(Symbol.(vec(header)))
    [NamedTuple{keys}(Tuple(row)) for row in eachrow(rows)]
end

function load_figures_data(folder)
    summary = Dict((r.case, r.method) => r for r in read_table(folder, "summary"))
    histories = Dict{Tuple{String,String},Vector{NamedTuple}}()
    for r in read_table(folder, "histories")
        push!(get!(Vector{NamedTuple}, histories, (r.case, r.method)), r)
    end
    poles = Dict{Tuple{String,String},Vector{NamedTuple}}()
    for r in read_table(folder, "verified_curve_poles")
        push!(get!(Vector{NamedTuple}, poles, (r.case, r.method)), r)
    end
    domains = Dict{String,Vector{ComplexF64}}()
    for r in read_table(folder, "domains")
        z = complex(parse(Float64, r.real), parse(Float64, r.imag))
        push!(get!(Vector{ComplexF64}, domains, r.case), z)
    end
    (;summary, histories, poles, domains)
end

number(r, key) = parse(Float64, getproperty(r, key))
is_returned(r) = r.returned == "true"
returned_index(history) = only(findall(is_returned, history))

function through_best_validation(history)
    index = findfirst(row -> row.retrospective_best_validation == "true", history)
    if index !== nothing && isfinite(number(history[index], :validation_error))
        return view(history, 1:index)
    end
    return history
end

method_label(method) = STYLES[method][1]

function convergence_legend!(position, methods; nbanks=1, crosses=false)
    elements = [[LineElement(color=STYLES[method][2], linestyle=STYLES[method][3], linewidth=1.3),
        MarkerElement(color=STYLES[method][2], marker=STYLES[method][4], markersize=4)]
        for method in methods]
    labels = method_label.(methods)
    if crosses
        push!(elements, [MarkerElement(color=:black, marker=:xcross, markersize=6)])
        push!(labels, "nonfinite error")
    end
    Legend(position, elements, labels; orientation=:horizontal, nbanks, framevisible=true,
        padding=(4, 4, 2, 2), labelsize=10)
end

function convergence_panel!(axis, histories, case, methods, xfield)
    displayed = [through_best_validation(histories[(case, method)]) for method in methods]
    positive = [number(row, :validation_error) for rows in displayed for row in rows
        if isfinite(number(row, :validation_error)) && number(row, :validation_error) > 0]
    lower = minimum(positive; init=1.0) / 10
    upper = maximum(positive; init=1.0) * 10
    has_nonfinite = false
    has_zero = false
    for (method, history) in zip(methods, displayed)
        _, color, linestyle, marker = STYLES[method]
        x = Int.(number.(history, xfield))
        errors = number.(history, :validation_error)
        y = [isfinite(e) ? max(e, lower) : NaN for e in errors]
        scatterlines!(axis, x, y; color, linestyle, marker, linewidth=1.0, markersize=2.8, alpha=0.65)
        bad = findall(!isfinite, errors)
        has_nonfinite |= !isempty(bad)
        has_zero |= any(iszero, errors)
        if !isempty(bad)
            scatter!(axis, x[bad], fill(upper, length(bad)); color, marker=:xcross, markersize=7)
        end
        if isfinite(last(errors))
            scatter!(axis, [last(x)], [last(y)];
                color=:transparent, strokecolor=:black, strokewidth=1.8, markersize=10)
        end
    end
    (; nonfinite=has_nonfinite, zero=has_zero)
end

function pole_coordinates(rows)
    z = complex.(number.(rows, :real), number.(rows, :imag))
    residues = number.(rows, :residue_abs)
    located = isfinite.(z) .& (getproperty.(rows, :finite_pole) .== "true")
    colored = located .& isfinite.(residues) .& (getproperty.(rows, :finite_residue) .== "true")
    (; z, residues, located, colored)
end

function pole_panel!(axis, rows, result; exact_degree=nothing)
    z, residues, located, colored = pole_coordinates(rows)
    inside = located .& (abs.(real.(z)) .<= POLE_WINDOW) .& (abs.(imag.(z)) .<= POLE_WINDOW)
    keep = inside .& colored .& (residues .> 0)
    zero_residue = inside .& colored .& iszero.(residues)
    unknown = inside .& .!colored
    any(keep) && scatter!(axis, real.(z[keep]), imag.(z[keep]); color=log10.(residues[keep]),
        colormap=Reverse(:thermal), colorrange=(-4, 0), markersize=5)
    any(zero_residue) && scatter!(axis, real.(z[zero_residue]), imag.(z[zero_residue]);
        color=:transparent, strokecolor=:black, strokewidth=0.8, markersize=5)
    any(unknown) && scatter!(axis, real.(z[unknown]), imag.(z[unknown]); color=:gray40,
        marker=:xcross, markersize=6)
    notes = String[]
    if result.stop_reason == "nonfinite_weight"
        if result.returned_block == "0" && result.failed_block == "1" && result.failed_leg == "2"
            push!(notes, "Endpoint pair failed\nNonfinite coefficient")
        else
            push!(notes, "Nonfinite coefficient\nblock $(result.failed_block), leg $(result.failed_leg)")
        end
    elseif result.extraction_status != "ok"
        push!(notes, "Candidate extraction unavailable")
    elseif !any(located)
        result.stop_reason == "max_degree" && push!(notes, "Degree cap")
        push!(notes, "No finite candidates")
    end
    outside = count(located .& .!inside)
    outside > 0 && push!(notes, "outside window: $outside")
    candidate_count = count(located)
    if exact_degree !== nothing && candidate_count != exact_degree
        push!(notes, "$candidate_count finite candidates, exact degree $exact_degree")
    end
    any(unknown) && push!(notes, "unavailable residues: $(count(unknown))")
    any(!, located) && push!(notes, "nonfinite candidates: $(count(!, located))")
    isempty(notes) || text!(axis, 0.03, 0.012; text=join(notes, "\n"), space=:relative,
        align=(:left, :bottom), fontsize=8)
    any(zero_residue) || any(unknown)
end
