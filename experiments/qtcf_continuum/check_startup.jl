# Check startup nodes and the exclusion of terminal greedy samples.
include("run.jl")

function certify_startup()
    fits = deserialize(joinpath(DEST,"results.jls"))
    cases = Dict(c.slug=>c for c in make_cases())
    starts = NamedTuple[]
    for ((slug, method), fit) in sort!(collect(fits);by=first)
        @assert method in ("QTCF", "forced_Classic", "native_squared")
        quadratic = method == "QTCF"
        c = cases[slug]
        squared = method == "native_squared"
        L = length(c.domain)
        s = quadratic ? L/2 : zero(Float64)
        first_support = point(c.domain,s)
        partner = c.reflection(first_support)
        terminal = point(c.domain,L)
        first_model = fit.best == 0 ? fit.partial : first(fit.history).interpolant
        initial_nodes = nodes(first_model)
        reflected_fraction = NaN
        if quadratic
            @assert initial_nodes == [first_support]
            @assert all(length(supports) == 2block-1 for (block,supports) in enumerate(fit.supports))
            # The first leg is greedy. The partner is forced by reflection.
            @assert all(supports[end-1] != terminal for supports in fit.supports[2:end])
        elseif squared
            @assert initial_nodes == [first_support^2]
            @assert partner == terminal
            @assert first(fit.supports)[1:2] == [first_support,partner]
            @assert first(fit.attempts).parameter == 0
            transformed_terminal = point(c.domain,L/2)^2
            @assert all(transformed_terminal ∉ nodes(rec.interpolant) for rec in fit.history)
            reflected_fraction = 1.0
        else
            @assert first(initial_nodes) == first_support
            @assert first(fit.attempts).parameter == 0
            if partner != first_support
                @assert partner == terminal
                @assert fit.attempts[2].parameter == L
                reflected_fraction = 1.0
                if fit.best == 0
                    @assert fit.reason === :nonfinite_weight && fit.failure.block == 1 &&
                        fit.failure.leg == 2 && fit.failure.node == partner
                    @assert length(initial_nodes) == 1
                else
                    @assert initial_nodes == [first_support,partner]
                end
            else
                @assert initial_nodes == [first_support]
            end
            n = length(initial_nodes)
            @assert all(length(supports) == n+2(block-1) for (block,supports) in enumerate(fit.supports))
            @assert all(0 < attempt.parameter < L for attempt in fit.attempts
                if attempt.block > 1 && attempt.leg == 1)
        end
        test_parameters,_ = collect(fit.path,:test)
        @assert all(0 < t < length(fit.path.path) for t in test_parameters)
        push!(starts,(;case=slug,method,parameter_fraction=s/L,
            reflected_parameter_fraction=reflected_fraction,
            initial_supports=fit.best == 0 ? 0 : length(first(fit.supports)),
            first_support_real=real(first_support),first_support_imag=imag(first_support),
            first_pair_completed=partner != first_support && fit.best > 0,
            terminal_greedy_sample=false))
    end
    table("startup_certificates",starts)
    println("Certified ", length(starts), " startup records")
end

abspath(PROGRAM_FILE)==(@__FILE__) && certify_startup()
