# Repeat both absolute-value fits with refinement 3/3 and certify the returned quotients.
# Run in a prepared frozen experiment, or supply its directory and an output directory.
const STRESS_EXPERIMENT = isempty(ARGS) ? (@__DIR__) : abspath(first(ARGS))
const STRESS_OUTPUT = length(ARGS) < 2 ? STRESS_EXPERIMENT : abspath(ARGS[2])
include(joinpath(STRESS_EXPERIMENT, "certify_intervals.jl"))
mkpath(STRESS_OUTPUT)
function fit_parameters(case, method::Symbol)
    options = merge(experiment_parameters(case.family, method),(;tol=REL_TOL))
    return method === :QTCF ? merge(options,(;initial_refinement=3,refinement=3)) : options
end
BLAS.set_num_threads(1)
global_logger(ConsoleLogger(stderr, Logging.Error))
rows, axis, bounds, brackets, comparisons = [NamedTuple[] for _ in 1:5]
for case in make_cases()
    case.slug in ("real","imaginary") || continue
    fit = production_fit(case;log_attempts=true)
    label = "QTCF"
    trace = filter(r -> r.case==case.slug && r.method==label, REDUCTIONS)
    reduction_peak = maximum(r.reduced_max for r in trace; init=0.0)
    target = validation_error(fit.r, fit.kind, case)
    raw = RFA.evaluate(fit.r,case.check)
    entry = (;r=fit.r,kind=:QTCF,physical_nodes=length(nodes(fit.r)))
    row=(;case=case.slug,method=label,fit.reason,error=target.error,nonfinite=target.nonfinite,
        raw_nonfinite=count(!isfinite,raw),nodes=entry.physical_nodes,blocks=length(fit.history),
        reduction_peak,failed_reductions=count(r -> r.outcome == "failed", trace))
    push!(comparisons,row)
    println("STRESS SUMMARY ",row)
    flush(stdout)
    certify_model(rows,axis,bounds,brackets,(case.slug,"refinement_3_3",label),entry,"")
    println("STRESS CERTIFICATE ",last(rows))
    flush(stdout)
    for (name,data) in (("stress_3_3",comparisons),("stress_3_3_certificates",rows),
        ("stress_3_3_axis",axis),("stress_3_3_bounds",bounds),("stress_3_3_brackets",brackets))
        save_cert(STRESS_OUTPUT,name*".csv",data)
    end
end
