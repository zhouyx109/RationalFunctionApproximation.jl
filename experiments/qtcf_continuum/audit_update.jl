# Check the normalized update on fixed nodes and compare it with exact arithmetic.
# Load certify_intervals.jl before this file. All fit and certification
# settings come from the frozen continuum experiment.

function replay_update(fit, method)
    source = fit.r
    values = source.vals
    z = nodes(source)
    r = QuadraticThiele(z[1:1], values[1:1]; block_sizes=[1], method)
    offset = 1
    for pair in source.pairs[2:end]
        count = pair[1] == pair[2] ? 1 : 2
        inds = offset .+ (1:count)
        RFA.add_node!(r, z[inds], values[inds], method)
        offset += count
    end
    @assert offset == length(z) && isequal(nodes(r), z)
    return r
end

# A short exact replay of the unscaled block formula, independent of the new
# representation. Raw exact interpolation can grow costly on long prefixes.
function canonical_prefix(source, blocks)
    pairs = map(pair -> exact.(pair), source.pairs[1:blocks])
    numerators = [pair[1] == pair[2] ? EQ[-pair[1],1] :
        pmul(EQ[-pair[1],1], EQ[-pair[2],1]) for pair in pairs]
    coefficients = [EQ[exact(first(source.vals))]]
    offset = 1
    for j in 2:blocks
        pair = pairs[j]
        count = pair[1] == pair[2] ? 1 : 2
        reduced = EQ[]
        for k in 1:count
            z = pair[k]
            # Use the projective pair over exact rationals. An intermediate
            # infinite quotient is valid when the final reduced value is finite.
            n, d = exact(source.vals[offset+k]), one(EQ)
            for l in 1:j-1
                n, d = peval(numerators[l], z)*d, n-peval(coefficients[l],z)*d
            end
            iszero(d) && error("canonical exact block $j has zero reduced denominator")
            push!(reduced, n/d)
        end
        if count == 1
            push!(coefficients, reduced)
        else
            slope = (last(reduced)-first(reduced))/(last(pair)-first(pair))
            push!(coefficients, EQ[first(reduced)-slope*first(pair),slope])
        end
        offset += count
    end
    return coefficients, numerators
end

function canonical_prefix_error(fit, case)
    blocks = min(6, fit.best)
    prefix = fit.history[blocks].interpolant
    bs, as = canonical_prefix(prefix, blocks)
    p, q, _ = exact_continuants(prefix)
    peak = 0.0
    checked = 0
    for index in round.(Int,range(2,length(case.check)-1;length=12))
        point = case.check[index]
        point in nodes(prefix) && continue
        z = exact(point)
        n, d = peval(last(bs), z), one(EQ)
        for j in blocks-1:-1:1
            n, d = peval(bs[j],z)*n+peval(as[j],z)*d, n
        end
        pn, qn = peval(p,z), peval(q,z)
        (iszero(d) || iszero(qn)) && continue
        expected, actual = n/d, pn/qn
        magnitude = max(one(QQ),abs2(expected),abs2(actual))
        peak = max(peak, sqrt(Float64(abs2(actual-expected)/magnitude)))
        checked += 1
    end
    return (;blocks, checked, peak)
end

function compare_update_grid(a, b, case)
    delta = 0.0
    bad = 0
    for z in case.check
        av, bv = lifted(a,z,:QTCF), lifted(b,z,:QTCF)
        if isfinite(av) && isfinite(bv)
            delta = max(delta, Float64(abs(av-bv)))
        else
            bad += 1
        end
    end
    return (;max_finite_difference=delta, points_with_nonfinite=bad)
end

function exact_same_quotient(a, b)
    # A differing cross product at one point disproves equality quickly.
    z = EQ(QQ(1,11),QQ(1,7))
    an, ad = exact_stored_pair(a,z)
    bn, bd = exact_stored_pair(b,z)
    an*bd == bn*ad || return false
    ap, aq, _ = exact_continuants(a)
    bp, bq, _ = exact_continuants(b)
    return ap == bp && aq == bq || pmul(ap,bq) == pmul(bp,aq)
end

function exact_stored_pair(r, z)
    r.chart !== nothing && return exact_chart_pair(r,z)
    bs, as = exact_blocks(r)
    n, d = peval(last(bs),z), one(EQ)
    for j in length(bs)-1:-1:1
        n, d = peval(bs[j],z)*n+peval(as[j],z)*d, n
    end
    return n, d
end

# Compare scalar evaluation with the exact stored quotient on a spaced set of
# independent probes. A removable 0/0 is recorded separately.
function stored_evaluation_gap(r, case)
    compared = 0
    indeterminate = 0
    failures = 0
    peak = 0.0
    for index in unique(round.(Int,range(1,length(case.check);length=12)))
        point = case.check[index]
        p, q = exact_stored_pair(r, exact(point))
        if iszero(q)
            indeterminate += Int(iszero(p))
            failures += Int(!iszero(p))
            continue
        end
        value = lifted(r,point,:QTCF)
        if isfinite(value)
            expected = p/q
            gap = abs2(exact(value)-expected)/max(one(QQ),abs2(expected))
            peak = max(peak,sqrt(Float64(gap)))
            compared += 1
        else
            failures += 1
        end
    end
    return (;compared,indeterminate,failures,peak)
end

function audit_update(fits, make_update=case -> NormalizedQTCF();
    original="QTCF", revised="normalized_QTCF", out=@__DIR__)
    mkpath(out)
    compare, failures = NamedTuple[], NamedTuple[]
    for case in make_cases()
        slug = case.slug
        before, after = fits[(slug,original)], fits[(slug,revised)]
        before.best > 0 && after.best > 0 || error("Missing fit for $slug")
        same_selected_nodes = isequal(nodes(before.r),nodes(after.r))
        adaptive_grid = compare_update_grid(before.r,after.r,case)
        fixed_nodes_completed = false
        fixed_exact = false
        fixed_grid = (;max_finite_difference=NaN,points_with_nonfinite=0)
        prefix = (;blocks=0,checked=0,peak=NaN)
        try
            fixed = replay_update(before, make_update(case))
            fixed_nodes_completed = true
            fixed_grid = compare_update_grid(before.r,fixed,case)
            fixed_exact = exact_same_quotient(before.r,fixed)
        catch e
            push!(failures,(;case=slug,method=revised,check="fixed_node_replay",error=sprint(showerror,e)))
        end
        try
            prefix = canonical_prefix_error(after,case)
        catch e
            push!(failures,(;case=slug,method=revised,check="exact_canonical_prefix",error=sprint(showerror,e)))
        end
        stored = stored_evaluation_gap(after.r,case)
        push!(compare,(;case=slug,original_nodes=length(nodes(before.r)),normalized_nodes=length(nodes(after.r)),
            same_selected_nodes,adaptive_max_finite_difference=adaptive_grid.max_finite_difference,
            adaptive_points_with_nonfinite=adaptive_grid.points_with_nonfinite,
            fixed_nodes_completed,fixed_exact_same_quotient=fixed_exact,
            fixed_max_finite_difference=fixed_grid.max_finite_difference,
            fixed_points_with_nonfinite=fixed_grid.points_with_nonfinite,
            exact_prefix_blocks=prefix.blocks,exact_prefix_points=prefix.checked,
            exact_prefix_max_scaled_difference=prefix.peak,
            stored_points_compared=stored.compared,stored_indeterminate_points=stored.indeterminate,
            stored_failed_points=stored.failures,stored_max_scaled_difference=stored.peak))
        println("UPDATE REPLAY ",last(compare))
        flush(stdout)
        save_cert(out,"update_replay.csv",compare)
        save_cert(out,"update_audit_failures.csv",failures)
    end
    return (;compare,failures)
end
