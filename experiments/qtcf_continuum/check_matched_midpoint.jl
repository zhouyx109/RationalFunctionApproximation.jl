# Compare the Group L fits with both methods starting at the same fixed node.
include("check_pole_quality.jl")

function matched_midpoint()
    BLAS.set_num_threads(1)
    global_logger(ConsoleLogger(stderr, Logging.Error))
    fits = deserialize(joinpath(DEST, "results.jls"))
    rows = NamedTuple[]
    for case in make_cases()
        case.family === :schwarz && case.parity === :none || continue
        midpoint = point(case.domain, length(case.domain) / 2)
        @assert midpoint == case.reflection(midpoint)
        scale = case_metrics(case).absolute_target / REL_TOL

        paired = paired_fit(case; midpoint_start=true, log_attempts=true)
        original = fits[(case.slug, "QTCF")]
        @assert first(nodes(paired.r)) == first(nodes(original.r)) == midpoint
        for (method, fit) in (("midpoint_paired_TCF", paired), ("QTCF", original))
            validation = validation_error(fit.r, fit.kind, case)
            denominator = exact_model(fit.r, fit.kind).q
            push!(rows, (;case=case.slug, method, nodes=length(nodes(fit.r)),
                reduced_denominator_degree=degree(denominator),
                relative_error=validation.error / scale, nonfinite=validation.nonfinite,
                stop_reason=fit.reason,
                exact_denominator_symmetry=denominator_has_exact_symmetry(denominator, case.reflection)))
            println(last(rows))
        end
    end
    table("matched_midpoint", rows)
end

abspath(PROGRAM_FILE) == (@__FILE__) && matched_midpoint()
