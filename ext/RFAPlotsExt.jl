module RFAPlotsExt

using RationalFunctionApproximation, ComplexRegions, Plots

# Plot a convergence curve, showing which steps had bad poles
function plot_convergence(r::RFA.AbstractApproximation)
    deg, err, _, allowed, best = get_history(r)
    fig = plot(xlabel="degree", ylabel="max error", yscale=:log10, legend=false)
    seriescolor = [ all(b) ? :darkblue : :red for b in allowed]
    scatter!(fig, [deg[best]], [err[best]],
        markercolor=RGBA(1,1,1,0), markerstrokecolor=:gold, msw=5, markersize=7 )
    scatter!(fig, deg, err; seriescolor, msw=0, markersize=5)
    return fig
end

function RFA.errorplot(r::RFA.AbstractApproximation; use_abs=false)
    if r isa RFA.DiscreteApproximation
        _, err = check(r; quiet=true)
        t = findall(r.test_index)
        xlabel = "sample index"
    else
        N = cld(1000, max(1, length(nodes(r))))
        t, _, err = check(r; refinement=N, quiet=true, prenodes=true)
        xlabel = "boundary parameter"
    end
    fig = plot(; xlabel, ylabel="error", legend=false)
    if use_abs
        plot!(fig, t, abs.(err))
        ylabel!(fig, "| error |")
    elseif isreal(err)
        plot!(fig, t, real.(err))
    else
        plot!(fig, t, real.(err))
        plot!(fig, t, imag.(err))
        ylabel!(fig, "Re, Im error")
    end
    return fig
end

function axisbox(z)
    xbox, ybox = ComplexRegions.enclosing_box(z, 1.42)
    center = [(box[2] + box[1]) / 2 for box in (xbox, ybox)]
    radius = maximum((box[2] - box[1]) / 2 for box in (xbox, ybox))
    xbox = center[1] .+ [-radius, radius]
    ybox = center[2] .+ [-radius, radius]
    return xbox, ybox
end

plot_poles(r::RFA.AbstractApproximation) = RFA.poleplot(r, 0)

function RFA.poleplot(r::RFA.AbstractApproximation, idx::Integer)
    fig = plot(xlabel="Re(z)", ylabel="Im(z)", aspect_ratio=1, legend=false)
    z = test_points(r; include_selected=true)
    if r isa RFA.DiscreteApproximation
        scatter!(fig, real(z), imag(z); markersize=2)
    else
        RFA.isclosed(domain(r)) && (z = [z; first(z)])
        plot!(fig, real(z), imag(z))
    end
    zp = iszero(idx) ? RFA.poles(r) : RFA.poles(rewind(r, idx))
    color = [RFA._pole_allowed(r.allowed, z) ? :black : :red for z in zp]
    scatter!(fig, real(zp), imag(zp); color, markersize=4)
    xbox, ybox = axisbox(z)
    xlims!(fig, xbox...)
    ylims!(fig, ybox...)
    return fig
end

end # module
