# Fit every target with the listed methods and save the timing, error, and weight records.
using RationalFunctionApproximation, ComplexRegions, Serialization, Statistics
using Logging, Printf, LinearAlgebra, SHA, Dates

const PACKAGE = joinpath(@__DIR__, "frozen", "RFA")
realpath(pathof(RFA)) == realpath(joinpath(PACKAGE,"src","RationalFunctionApproximation.jl")) ||
    error("The local RFA package is required")
# Read the geometry and settings from the archived experiment.
parameters=read(joinpath(PACKAGE,"test","qtcf_parameters.jl"),String)
include_string(@__MODULE__,first(split(parameters,"# Figure labels and display.")),"paper_parameters")
experiments=read(joinpath(PACKAGE,"test","qtcf_test.jl"),String)
geometry=first(split(last(split(experiments,"# Interval data and independent validation points.")),
    "include(\"qtcf_output.jl\")"))
include_string(@__MODULE__,geometry,"paper_geometry")
set_eval_method(EXPERIMENT_EVAL_METHOD)
set_weight_method(EXPERIMENT_WEIGHT_METHODS.TCF)

function run_source_hashes(root, project)
    files=String[]
    for folder in ("src","ext","test")
        for (directory,subdirs,names) in walkdir(joinpath(root,folder))
            filter!(s -> s ∉ ("qtcf_output","qtcf_discrete_output"),subdirs)
            append!(files,(joinpath(directory,name) for name in names
                if endswith(name,".jl") || name in ("Project.toml","Manifest.toml")))
        end
    end
    append!(files,(project,joinpath(root,"Project.toml")))
    Dict(relpath(file,root)=>bytes2hex(sha256(read(file))) for file in unique(files))
end

function curve_points(domain,n;periodic=false)
    parameters=range(0,length(domain);length=n+periodic)
    periodic && (parameters=parameters[1:end-1])
    [point(domain,s) for s in parameters]
end

function validation_grid(case)
    (;points=curve_points(case.fit_domain,SCHWARZ_AUDIT_N;
        periodic=RFA.isclosed(case.fit_domain)))
end

include("traced_classic.jl")

const DEST = @__DIR__
const REL_TOL = 4e-11
const REDUCTIONS = NamedTuple[]
const EVALUATOR_CHECKS = NamedTuple[]
const SUMMARY = NamedTuple[]
const HISTORIES = NamedTuple[]
const POLE_ROWS = NamedTuple[]
const DIAGNOSTICS = NamedTuple[]
const DOMAIN_ROWS = NamedTuple[]
const SUPPORT_ROWS = NamedTuple[]
const RESULTS = Dict{Tuple{String,String},Any}()
const PHYSICAL_CAP = 399
const DENOMINATOR_CAP = 200

# Save independent experiment tables without replacing the paper outputs.
function table(name, rows)
    isempty(rows) && return
    fields = propertynames(first(rows))
    open(joinpath(DEST, name * ".csv"), "w") do io
        println(io, join(fields, ','))
        for row in rows
            print(io, join((begin
                value = string(getproperty(row, field))
                occursin(r"[,\"\n]", value) ? "\"" * replace(value, "\"" => "\"\"") * "\"" : value
            end for field in fields), ','), '\n')
        end
    end
end

function save_tables()
    for (name, rows) in (("summary", SUMMARY), ("histories", HISTORIES),
        ("poles", POLE_ROWS), ("diagnostics", DIAGNOSTICS),
        ("domains", DOMAIN_ROWS), ("supports", SUPPORT_ROWS))
        table(name, rows)
    end
end

# All errors are measured in the original fitting coordinate.
function lifted(r, z, kind)
    kind in (:forced, :QTCF) && return RFA.evaluate(r, z, OneDiv())
    isempty(r.weights) && return zero(z)
    value = RFA.evaluate(r, z*z, OneDiv())
    return kind === :even ? value : z * value
end

function validation_error(r, kind, case)
    error = 0.0
    nonfinite = 0
    for (z, target) in zip(case.check, case.f.(case.check))
        value = abs(target - lifted(r, z, kind))
        if !isfinite(value)
            nonfinite += 1
            error = Inf
        else
            error = max(error, value)
        end
    end
    (; error, nonfinite)
end

function case_metrics(case)
    check = case.check
    diameter = if case.family === :interval
        2.0
    else
        maximum(abs(z-w) for z in check for w in check)
    end
    (; diameter, absolute_target=REL_TOL*case.target_scale)
end

function make_cases()
    cases = Any[]
    for c in interval_cases
        push!(cases, (; slug=c.slug, title=c.name, f=c.f, domain=c.domain,
            reflection=z -> -z, value_map=(z,y,w) -> y,
            rotation=1.0+0.0im, parity=:even, family=:interval,
            check=c.points, parameters=EXPERIMENT_PARAMETERS.interval))
    end
    for c in qtcf_cases
        mapvalue = c.function_reflection === (-) ? ((z,y,w) -> -y) :
            ((z,y,w) -> c.function_reflection(y))
        push!(cases, (; slug=c.slug, title=c.name, f=conj, domain=c.fit_domain,
            reflection=c.reflection, value_map=mapvalue, rotation=c.rotation,
            parity=c.reflection === (-) ? :odd : :none, family=:schwarz,
            check=validation_grid(c).points, parameters=EXPERIMENT_PARAMETERS.schwarz))
    end
    push!(cases, (; slug="odd_sine", title="sin(x) on [-1,1]", f=sin,
        domain=Segment{Float64}(-1, 1), reflection=z -> -z,
        value_map=(z,y,w) -> -y, rotation=1.0+0.0im, parity=:odd,
        family=:interval, check=abs_check, parameters=EXPERIMENT_PARAMETERS.interval))
    map(cases) do case
        center = point(case.domain, length(case.domain)/2)
        check = case.parity === :none || any(iszero, case.check) ? case.check : vcat(case.check, [center])
        merge(case, (;check, target_scale=maximum(abs, case.f.(check))))
    end
end

function fit_parameters(case, method)
    merge(experiment_parameters(case.family, method), (;tol=REL_TOL))
end

# Count original target calls separately from represented physical supports.
mutable struct CountedTarget{F} <: Function
    f::F
    calls::Int
    unique_coordinates::Set{Any}
end
CountedTarget(f) = CountedTarget(f,0,Set{Any}())
function (f::CountedTarget)(z)
    f.calls += 1
    canonical = z isa Complex ? complex(iszero(real(z)) ? zero(real(z)) : real(z),
        iszero(imag(z)) ? zero(imag(z)) : imag(z)) : (iszero(z) ? zero(z) : z)
    push!(f.unique_coordinates, canonical)
    f.f(z)
end
counter_snapshot(f::CountedTarget) = (;target_calls=f.calls, unique_target_calls=length(f.unique_coordinates))
counter_snapshot(f::Function) = (;target_calls=-1, unique_target_calls=-1)
classify_weight(w) = isnan(w) ? :NaN : isfinite(w) ? :finite : :Inf
function attempt_row(case,method,kind,block,leg,r,z,y,w;accepted=false,reason="",parameter=NaN)
    magnitude=abs(w)
    finite=isfinite(w)
    before=isempty(nodes(r)) ? NaN : abs(y-RFA.evaluate(r,z,OneDiv()))
    (;case=case.slug,method,kind=string(kind),attempt=length(nodes(r))+1,block,leg,
      prefix_nodes=length(nodes(r)),parameter,node_real=real(z),node_imag=imag(z),
      target_real=real(y),target_imag=imag(y),weight_real=real(w),weight_imag=imag(w),
      weight_abs=magnitude,log10_weight_abs=iszero(magnitude) ? -Inf : log10(magnitude),
      classification=string(classify_weight(w)),large_1e8=finite&&magnitude>=1e8,
      large_1e12=finite&&magnitude>=1e12,large_1e16=finite&&magnitude>=1e16,
      pre_error=before,accepted,reason=string(reason),native_guard_would_accept=!isnan(w))
end

# Choose a point by its own error, then insert it and its reflection one at a time.
# The main baseline starts at an endpoint. The separate control starts at the midpoint.
function paired_fit(case;log_attempts=false,method_label="forced_Classic",midpoint_start::Bool=false)
    kind = :forced
    coefficient = EXPERIMENT_WEIGHT_METHODS.TCF
    d,f,reflection=case.domain,(log_attempts ? CountedTarget(case.f) : case.f),case.reflection
    options=fit_parameters(case,:QTCF)
    refinement=options.refinement;initial_refinement=options.initial_refinement
    endpoint_start = !midpoint_start
    L=Float64(length(d));start=endpoint_start ? 0.0 : L/2
    z0=point(d,start);y0=f(z0);S=promote_type(typeof(z0),typeof(y0))
    mate=reflection(z0)
    startup_pair=mate!=z0
    max_iter=options.max_iter-Int(startup_pair)
    sampling_rule=endpoint_start ? "first endpoint and its reflection; interior greedy mesh" :
        "midpoint startup; original-variable residual and interior block mesh"
    r=Thiele(S[z0],S[y0],S[y0])
    attempts=NamedTuple[]
    if log_attempts
        empty=Thiele(S[],S[],S[])
        push!(attempts,attempt_row(case,method_label,kind,1,1,empty,z0,y0,y0;
            accepted=true,reason="startup",parameter=start))
    end
    path=RFA.DiscretizedPath(d,Float64[0,L];refinement=max(refinement,initial_refinement),
        maxpoints=max(2max_iter+2,max_iter*refinement))
    selected_params=Float64[start];selected_nodes=S[z0];selected_values=S[y0]
    if startup_pair
        s=mate==point(d,L) ? L : RFA._reflection_param(d,start,reflection)
        raw=f(mate);v=case.value_map(z0,y0,mate)
        weight=RFA._new_weight(coefficient,r.nodes,r.weights,mate,v)
        finite=isfinite(weight)
        log_attempts && push!(attempts,attempt_row(case,method_label,kind,1,2,r,mate,v,weight;
            accepted=finite,reason=finite ? "startup" : "nonfinite_weight",parameter=s))
        if !finite
            # The initial pair is atomic. A single support is a diagnostic partial model.
            partial=copy(r)
            empty=Thiele(S[],S[],S[])
            failure=(;block=1,leg=2,node=mate,weight)
            return (;r=empty,history=RFA.IterationRecord[],best=0,reason=:nonfinite_weight,
                supports=Vector{S}[],kind,path,production=nothing,counts=counter_snapshot(f),
                attempts,failure,partial,residual_evaluations=0,center_limit_supplied=false,
                sampling_rule)
        end
        push!(r.nodes,mate);push!(r.values,v);push!(r.weights,weight)
        push!(selected_params,s);push!(selected_nodes,mate);push!(selected_values,raw)
    end
    check_center=case.parity!==:none
    center=check_center ? point(d,L/2) : z0
    center_value=check_center ? f(center) : y0
    RFA.reset!(path,RFA._mesh_params(d,(),selected_params,Float64);refinement=path.nref)
    τ,num_ref=path.points,path.nref;rows=size(τ,1)
    active=Int[]
    fτ=Matrix{S}(undef,size(τ))
    fmax=0.0
    history=RFA.IterationRecord[RFA.IterationRecord(r,NaN,missing)]
    supports=[copy(selected_nodes)];errors=Float64[]
    ztest,ytest=S[],S[];out,a,b=S[],S[],S[]
    reason,best=:max_degree,1;failure=nothing;partial=nothing;residual_evaluations=0
    while true
        empty!(active)
        for j in 2:num_ref+1, i in 1:length(path.next)
            index=(j-1)*rows+i
            0 < path.params[index] < L && τ[index] ∉ selected_nodes && push!(active,index)
        end
        RFA._sweep!(fτ,f,τ,active)
        empty!(ztest);empty!(ytest)
        append!(ztest,view(τ,active));append!(ytest,view(fτ,active))
        for i in eachindex(selected_params)
            selected_params[i] < L || continue
            push!(ztest,selected_nodes[i]);push!(ytest,selected_values[i])
        end
        if check_center
            push!(ztest,center)
            push!(ytest,center_value)
        end
        for buffer in (out,a,b);resize!(buffer,length(ztest));end
        resize!(errors,length(ztest))
        RFA.evaluate!(out,r,ztest,a,b)
        for i in eachindex(out)
            isnan(out[i]) && (out[i] = lifted(r,ztest[i],kind))
        end
        @. errors=abs(ytest-out)
        residual_evaluations+=length(ztest)
        history[end].error=maximum(errors)
        fmax=max(fmax,maximum(abs,ytest))
        stopping_scale=case.parity===:none ? fmax : case.target_scale
        reason,best=RFA.quitting_check(history,options.stagnation,options.tol,stopping_scale,max_iter,true)
        reason===:iterating || break
        eligible=[k for k in eachindex(active) if τ[active[k]] ∉ selected_nodes && τ[active[k]]!=reflection(τ[active[k]])]
        if isempty(eligible);reason,best=:exhausted,RFA.best_acceptable(history,true);break;end
        _,worst=findmax(view(errors,eligible));index=active[eligible[worst]]
        sp,z,y=path.params[index],τ[index],fτ[index]
        w=reflection(z);partner=RFA._reflection_param(d,sp,_ -> w;candidate=L-sp)
        vraw=f(w);v=case.value_map(z,y,w)
        if !all(isfinite,(y,vraw,v))
            reason,best=:nonfinite_target,RFA.best_acceptable(history,true)
            failure=(;block=length(history)+1,leg=0,node=w,weight=oftype(v,NaN));break
        end
        trial=copy(r)
        for (leg,node,value,param) in ((1,z,y,sp),(2,w,v,partner))
            weight=RFA._new_weight(coefficient,trial.nodes,trial.weights,node,value)
            finite=isfinite(weight)
            log_attempts && push!(attempts,attempt_row(case,method_label,kind,length(history)+1,leg,
                trial,node,value,weight;accepted=finite,reason=finite ? "accepted" : "nonfinite_weight",parameter=param))
            if !finite
                failure=(;block=length(history)+1,leg,node,weight)
                partial=copy(trial)
                break
            end
            # Use the already-computed native coefficient; no duplicate weight calculation.
            push!(trial.nodes,node);push!(trial.values,value);push!(trial.weights,weight)
        end
        if failure!==nothing
            reason,best=:nonfinite_weight,RFA.best_acceptable(history,true);break
        end
        r=trial
        push!(history,RFA.IterationRecord(r,NaN,missing))
        append!(selected_params,(sp,partner));append!(selected_nodes,(z,w));append!(selected_values,(y,vraw))
        push!(supports,copy(selected_nodes))
        try
            num_ref=max(refinement,num_ref-1)
            RFA.reset!(path,RFA._mesh_params(d,(),selected_params,Float64);refinement=num_ref)
        catch exception
            exception isa Union{ArgumentError,KeyError,RFA.MaxRefinementException} || rethrow()
            reason,best=:refinement,RFA.best_acceptable(history,true);break
        end
    end
    best=RFA.best_acceptable(history,true)
    (;r=copy(history[best].interpolant),history,best,reason,supports,kind,path,
      production=nothing,counts=counter_snapshot(f),attempts,failure,partial,residual_evaluations,
      center_limit_supplied=false,sampling_rule)
end

# An external delegating method traces native scalar attempts in a separate untimed run.
# Timed native fits use uninstrumented Classic construction.
struct LoggedWeight{M,C} <: RFA.ThieleMethod
    method::M
    case::C
    label::String
    kind::Symbol
    attempts::Vector{NamedTuple}
end
function RFA._new_weight(m::LoggedWeight,z,w,node,value)
    coefficient=RFA._new_weight(m.method,z,w,node,value)
    prefix=Thiele(copy(z),fill(zero(eltype(w)),length(z)),copy(w))
    row=attempt_row(m.case,m.label,m.kind,length(z)+1,1,prefix,node,value,coefficient;
        accepted=!isnan(coefficient),reason=isnan(coefficient) ? "native_NaN_rejection" : "native_acceptance")
    # The native weight hook receives nodes and weights, not prefix target values.
    push!(m.attempts,merge(row,(;pre_error=NaN)))
    coefficient
end
function native_call(call,case,label,kind,log_attempts)
    attempts=NamedTuple[]
    if !log_attempts
        return call(),attempts
    end
    prior=RFA.default_weight_method()
    trace=LoggedWeight(prior,case,label,kind,attempts)
    try
        RFA.set_weight_method(trace)
        result=Base.invokelatest(call)
        r=first(result)
        if !isempty(nodes(r))
            empty=Thiele(eltype(r)[],eltype(r)[],eltype(r)[])
            pushfirst!(attempts,attempt_row(case,label,kind,1,1,empty,first(nodes(r)),first(values(r)),first(weights(r));
                accepted=true,reason="native_startup",parameter=0.0))
        end
        return result,attempts
    finally
        RFA.set_weight_method(prior)
    end
end

# Run the original midpoint-start QTCF with the package's Classic construction.
function production_fit(case;log_attempts=false)
    f = log_attempts ? CountedTarget(case.f) : case.f
    trace = log_attempts ? NamedTuple[] : nothing
    options = fit_parameters(case,:QTCF)
    startup = point(case.domain,length(case.domain)/2)
    weight_method = trace === nothing ? Classic() : TracedClassic(trace)
    comparison_scale = case.parity === :none ? nothing : case.target_scale
    result = approximate(f,case.domain,QTCF();options...,weight_method,eval_method=OneDiv(),
        reflection=case.reflection,value_reflection=case.value_map,
        initial_node=startup,pairs_only=true,comparison_scale,
        parameter_reflection=s -> length(case.domain)-s)
    if trace !== nothing
        append!(REDUCTIONS, [(;case=case.slug,method="QTCF",row...) for row in trace])
    end
    (;r=get_function(result),history=result.history,best=status(result).best,reason=status(result).reason,
      supports=[copy(nodes(h.interpolant)) for h in result.history],
      kind=:QTCF,path=result.path,production=result,
      counts=counter_snapshot(f),attempts=NamedTuple[],failure=nothing,partial=nothing,
      residual_evaluations=-1,center_limit_supplied=false,
      sampling_rule="native QTCF midpoint and complete pairs")
end

# Native transformed fits retain native TCF mesh, transformed residual and stopping.
function native_square_fit(case;log_attempts=false,method_label="native_squared")
    d=case.domain;f=log_attempts ? CountedTarget(case.f) : case.f;L=length(d)
    # Follow the original left half from its first endpoint to the midpoint.
    h=s -> point(d,L*s/2)
    left=h(0.0)
    u0=left^2
    representative(u)=u==u0 ? complex(left) : sqrt(complex(u))
    curve=Curve(s -> h(s)^2,s -> L*h(s)*ComplexRegions.tangent(d,L*s/2),0,1)
    g=if case.parity===:even
        u -> f(representative(u))
    elseif case.slug=="odd_sine"
        u -> iszero(u) ? 1.0+0im : f(representative(u))/representative(u)
    else
        tangent=ComplexRegions.tangent(d,L/2)
        u -> iszero(u) ? conj(tangent)/tangent : f(representative(u))/representative(u)
    end
    options=merge(fit_parameters(case,:TCF),(;max_iter=PHYSICAL_CAP÷2))
    call=() -> begin
        result=approximate(g,curve,TCF();options...)
        (get_function(result),result)
    end
    (r,result),attempts=native_call(call,case,method_label,case.parity,log_attempts)
    supports=map(result.history) do rec
        represented=ComplexF64[]
        for t in nodes(rec.interpolant)
            z=representative(t)
            push!(represented,z,-z)
        end
        case.parity===:odd && push!(represented,zero(ComplexF64))
        unique(represented)
    end
    (;r,history=result.history,best=status(result).best,reason=status(result).reason,
      supports,kind=case.parity,path=result.path,production=result,counts=counter_snapshot(f),
      attempts,failure=nothing,partial=nothing,residual_evaluations=-1,
      center_limit_supplied=case.parity===:odd,
      sampling_rule="squared image of the first endpoint; native transformed residual and interior mesh")
end

function original_degrees(r, kind)
    if isempty(nodes(r))
        return (0,0)
    end
    p,q=degrees(r)
    kind === :even && return (2p,2q)
    kind === :odd && return (2p+1,2q)
    return (p,q)
end

# Lift every finite scalar pencil candidate, retaining central and cancelled candidates.
function extracted(r,kind)
    isempty(nodes(r)) && return ComplexF64[],ComplexF64[]
    if kind in (:forced,:QTCF)
        return residues(r)
    end
    root, res = residues(r)
    pole, result = ComplexF64[],ComplexF64[]
    for (λ,A) in zip(root,res)
        if iszero(λ)
            push!(pole,0)
            push!(result,kind === :even ? 0 : A)
        else
            z=sqrt(complex(λ))
            append!(pole,(z,-z))
            append!(result,kind === :even ? (A/(2z),-A/(2z)) : (A/2,A/2))
        end
    end
    pole,result
end

function max_symmetry_defect(pole,reflection)
    finite=filter(isfinite,pole)
    isempty(finite) && return NaN
    maximum(minimum(abs(reflection(z)-w) for w in finite) for z in finite)
end

# Check that timing and tracing produce the same stored block coefficients.
function same_interpolant(a, b)
    isequal(nodes(a), nodes(b)) && isequal(weights(a), weights(b)) &&
        (!(a isa QuadraticThiele) || isequal(a.apolys, b.apolys))
end

function same_fit(a, b)
    same_interpolant(a.r, b.r) && a.best == b.best && a.reason == b.reason &&
        length(a.history) == length(b.history) &&
        all(isequal(x.error, y.error) && same_interpolant(x.interpolant, y.interpolant)
            for (x, y) in zip(a.history, b.history))
end
function audit_evaluators(case, r, method)
    z = case.check
    targets = case.f.(z)
    raw = RFA.evaluate(r,z)
    scalar = lifted.(Ref(r),z,Ref(:QTCF))
    for (evaluator,v) in (("stock_array",raw),("stock_scalar",scalar))
        finite = isfinite.(v)
        comparable = finite .& isfinite.(scalar)
        push!(EVALUATOR_CHECKS,(;case=case.slug,method,evaluator,points=length(z),nonfinite=count(!,finite),
            max_error=all(finite) ? maximum(abs,v.-targets) : Inf,
            max_difference_from_scalar=any(comparable) ? maximum(abs,v[comparable].-scalar[comparable]) : NaN,
            false_zeros=count(i -> iszero(v[i]) && isfinite(scalar[i]) && !iszero(scalar[i]),eachindex(v))))
    end
end

# Warm up before timing. Log weights and count target calls in a separate run.
function measured_fit(case,case_scale,method,runner)
    print(case.slug," ",method," warmup ... ");flush(stdout)
    warm=runner(false);times=Float64[];result=warm
    for _ in 1:3
        GC.gc();started=time_ns();result=runner(false);push!(times,(time_ns()-started)/1e9)
        same_fit(warm,result) || error("Nondeterministic fit $(case.slug)/$method")
    end
    traced=runner(true)
    same_fit(result,traced) || error("Attempt instrumentation changed fit $(case.slug)/$method")
    result=traced
    append!(DIAGNOSTICS,result.attempts)
    returned=result.best>0
    checker = validation_error
    check=returned ? checker(result.r,result.kind,case) : (;error=NaN,nonfinite=0)
    returned && result.kind===:QTCF && audit_evaluators(case,result.r,method)
    p,q=returned ? original_degrees(result.r,result.kind) : (-1,-1)
    pole,res=ComplexF64[],ComplexF64[];extraction_times=Float64[]
    extract_status=returned ? "ok" : "no_returned_fit"
    if returned
        try
            extracted(result.r,result.kind)
            for _ in 1:3
                started=time_ns();pole,res=extracted(result.r,result.kind);push!(extraction_times,(time_ns()-started)/1e9)
            end
        catch exception
            extract_status=string(typeof(exception))*": "*sprint(showerror,exception)
        end
    end
    physical=ComplexF64.(pole .* case.rotation);physical_res=ComplexF64.(res .* abs2(case.rotation))
    reflection=z -> case.rotation*case.reflection(z/case.rotation)
    defect=max_symmetry_defect(physical,reflection)
    support=returned ? unique(result.supports[result.best]) : nodes(result.r)
    hist_checks=[checker(rec.interpolant,result.kind,case) for rec in result.history]
    hist_errors=getproperty.(hist_checks,:error)
    best_error,best_validation=returned ? findmin(hist_errors) : (NaN,0)
    first_hit=findfirst(<=(case_scale.absolute_target),hist_errors)
    first_hit_nodes=first_hit===nothing ? 0 : length(unique(result.supports[first_hit]))
    first_hit_degree=first_hit===nothing ? -1 : last(original_degrees(result.history[first_hit].interpolant,result.kind))
    center=zero(eltype(case.check))
    center_error=returned ? abs(case.f(center)-lifted(result.r,center,result.kind)) : NaN
    row=(;case=case.slug,method,error=abs(case.rotation)*check.error,fit_seconds=median(times),
        fit_seconds_min=minimum(times),fit_seconds_max=maximum(times),
        extraction_seconds=isempty(extraction_times) ? NaN : median(extraction_times),extraction_status=extract_status,
        p_degree=p,q_degree=q,nodes=length(support),transformed_supports=length(nodes(result.r)),
        blocks=length(result.history),stop_reason=result.reason,returned_block=result.best,
        best_validation_error=abs(case.rotation)*best_error,best_validation_block=best_validation,
        pole_count=length(pole),symmetry_defect=defect,symmetry_defect_scaled=defect/case_scale.diameter,
        validation_diameter=case_scale.diameter,validation_target=case_scale.absolute_target,
        first_validation_hit_block=something(first_hit,0),first_validation_hit_nodes=first_hit_nodes,
        first_validation_hit_q_degree=first_hit_degree,returned_center_error=center_error,
        max_history_nodes=maximum(length, result.supports;init=0),
        nonfinite_evaluations=check.nonfinite,nonfinite_poles=count(!isfinite,physical),
        nonfinite_residues=count(!isfinite,physical_res),target_calls=result.counts.target_calls,
        unique_target_calls=result.counts.unique_target_calls,residual_evaluations=result.residual_evaluations,
        center_limit_supplied=result.center_limit_supplied,sampling_rule=result.sampling_rule,
        failed_block=result.failure===nothing ? 0 : result.failure.block,
        failed_leg=result.failure===nothing ? 0 : result.failure.leg,
        partial_error=result.partial===nothing ? NaN : validation_error(result.partial,result.kind,case).error)
    push!(SUMMARY,row)
    for (i,(z,A)) in enumerate(zip(physical,physical_res))
        push!(POLE_ROWS,(;case=case.slug,method,index=i,real=real(z),imag=imag(z),
            residue_real=real(A),residue_imag=imag(A),residue_abs=abs(A),finite_pole=isfinite(z),finite_residue=isfinite(A),
            lifted_central_multiple_pole=result.kind===:even&&iszero(z),
            simple_residue_does_not_describe_full_principal_part=result.kind===:even&&iszero(z)))
    end
    for (j,rec) in enumerate(result.history)
        p,q=original_degrees(rec.interpolant,result.kind)
        n=length(unique(result.supports[j]))
        @assert n<=PHYSICAL_CAP && q<=DENOMINATOR_CAP
        push!(HISTORIES,(;case=case.slug,method,block=j,nodes=n,
            transformed_supports=length(nodes(rec.interpolant)),p_degree=p,q_degree=q,
            sampled_error=rec.error,validation_error=abs(case.rotation)*hist_errors[j],returned=j==result.best,
            retrospective_best_validation=j==best_validation,nonfinite_evaluations=hist_checks[j].nonfinite,
            meets_validation_target=hist_errors[j]<=case_scale.absolute_target))
    end
    for (j,z) in enumerate(support)
        value=case.f(z)
        push!(SUPPORT_ROWS,(;case=case.slug,method,index=j,real=real(z),imag=imag(z),value_real=real(value),value_imag=imag(value)))
    end
    RESULTS[(case.slug,method)]=result
    save_tables();table("reductions",REDUCTIONS);table("evaluator_checks",EVALUATOR_CHECKS)
    serialize(joinpath(DEST,"results.jls"),RESULTS)
    @printf("error %.3e, degree (%d,%d), nodes %d, %.4f s, %s\n",row.error,row.p_degree,row.q_degree,row.nodes,row.fit_seconds,string(row.stop_reason))
    flush(stdout);row
end

function main()
    BLAS.set_num_threads(1);global_logger(ConsoleLogger(stderr,Logging.Error))
    source_before=run_source_hashes(PACKAGE,Base.active_project())
    open(joinpath(DEST,"provenance.txt"),"w") do io
        println(io,"Started UTC: ",Dates.now(Dates.UTC))
        println(io,"Frozen package: ",PACKAGE,"; Julia ",VERSION,"; BLAS threads 1; Float64")
        println(io,"QTCF and Paired TCF: ",INITIAL_REFINEMENT,"/",STEADY_REFINEMENT," once per complete post-startup block, with full mesh rebuilds.")
        println(io,"Native square TCF: ",SCALAR_INITIAL_REFINEMENT,"/",STEADY_REFINEMENT," once per inserted scalar support.")
        println(io,"Float64 relative tolerance ",REL_TOL," on all cases and for all methods. Group L block-based fits use the historical running maximum sampled target scale. Group P QTCF and Paired TCF use a common physical validation scale.")
        println(io,"All methods use Classic construction and OneDiv evaluation.")
        println(io,"Source base e19e600b4d0bca9ed12d6e6fa6369923ee6329b6. The experiment patch supplies a common Group P stopping scale for QTCF.")
        println(io,"Paired TCF starts at the first endpoint and inserts its reflection immediately unless fixed. A failed initial pair returns no model.")
        println(io,"Native square TCF starts at the square of the first endpoint. Its terminal parameter is the original midpoint.")
        println(io,"Common nominal caps: ",PHYSICAL_CAP," represented physical supports and denominator degree ",DENOMINATOR_CAP)
        println(io,"Reported relative errors use the maximum target magnitude on independent physical checks.")
        println(io,"Greedy sweeps use interior mesh points. Zero is checked during Group P Paired TCF fitting and on every Group P validation grid.")
        println(io,"The independent validation adds exact zero to each S-curve grid, yielding 1801 points.")
        println(io,"Full fits: one warmup, three measured repetitions, one untimed attempt trace; all histories checked deterministic")
        println(io,"Timing runs use original uninstrumented target functions. Counts from separate deterministic traced fit exclude independent validation; signed zeros canonicalized")
        println(io,"No rejection by coefficient magnitude; finite coefficients accepted; failures retain partial models")
        println(io,"Native scalar pole and residue candidates are retained without magnitude or location filters.")
        for key in sort!(collect(keys(source_before)));println(io,key," ",source_before[key]);end
    end
    requested = Set(ARGS)
    for case in make_cases()
        isempty(requested) || case.slug in requested || continue
        case_scale=case_metrics(case)
        measure(method,runner) = measured_fit(case,case_scale,method,runner)
        for (i,z) in enumerate(curve_points(case.domain,700))
            value=z*case.rotation;push!(DOMAIN_ROWS,(;case=case.slug,index=i,real=real(value),imag=imag(value)))
        end
        measure("QTCF",log -> production_fit(case;log_attempts=log))
        measure("forced_Classic",log -> paired_fit(case;log_attempts=log))
        if case.parity !== :none
            measure("native_squared",log -> native_square_fit(case;log_attempts=log))
        end
    end
    source_before==run_source_hashes(PACKAGE,Base.active_project()) || error("Frozen source changed")
    println("COMPLETED ",length(SUMMARY)," fits; frozen sources unchanged")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
