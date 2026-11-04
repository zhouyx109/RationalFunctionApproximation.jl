struct MaxRefinementException <: Exception
    msg::String
end

struct NaNException <: Exception
    msg::String
end

mutable struct DiscretizedPath{T,F,P}
    path::T              # original Curve or Path
    points::Matrix{P}    # first column = active points, others are refinements
    params::Matrix{F}    # first column = active points, others are refinements
    nref::Int            # number of refinements between nodes
    next::Vector{Int}    # next[i] = index of the next point in the path
end
"""
    DiscretizedPath(path, s::AbstractVector; kwargs...)
    DiscretizedPath(path, n::Integer=0; kwargs...)

Discretize a path, keeping the option of future making local refinements.

# Arguments
- `path`: a ComplexCurve or ComplexPath
- `s`: a vector of parameter values for the nodes
- `n`: number of equally spaced nodes on the path

# Keyword arguments
- `refinement`: number of refinements to make between consecutive points
- `maxpoints`: maximum number of points ever allowed

See also [`collect`](@ref), [`add_node!`](@ref).
"""
function DiscretizedPath(path::ComplexCurveOrPath, s::AbstractVector; refinement=0, maxpoints=length(s))
    n = length(s)
    n == 1 && throw(ArgumentError("Supply both endpoints or an empty parameter vector"))
    maxpoints >= n-1 || throw(ArgumentError("Number of nodes exceeds allocated size"))
    F = real_type(path)
    params = Matrix{F}(undef, maxpoints, refinement+1)
    for i in 1:n-1
        δ = (s[i+1] - s[i]) / (refinement + 1)
        for j in 0:refinement
            params[i, j+1] = s[i] + j * δ
        end
    end
    points = Matrix{typeof(point(path, 0))}(undef, maxpoints, refinement+1)
    for i in 1:n-1, j in 1:refinement+1
        points[i, j] = path(params[i, j])
    end
    next = n == 0 ? Int[] : [collect(2:n-1); 0]
    sizehint!(next, maxpoints)
    return DiscretizedPath{typeof(path),F,eltype(points)}(path, points, params, refinement, next)
end

function DiscretizedPath(path::ComplexCurveOrPath=ComplexRegions.Circle(0,1), n::Integer=0; kwargs...)
    DiscretizedPath(path, range(0, length(path), n); kwargs...)
end

function DiscretizedPath(p::DiscretizedPath, refinement::Int)
    s, _ = collect(p, :nodes)
    maxpoints = max(length(s), size(p.points, 1))
    return DiscretizedPath(p.path, s; refinement, maxpoints)
end

function reset!(d::DiscretizedPath, s::AbstractVector; refinement=0)
    if refinement + 1 > size(d.params, 2)
        throw(ArgumentError("New refinement level exceeds allocated size"))
    elseif length(s) > size(d.params, 1)
        throw(ArgumentError("New number of nodes exceeds allocated size"))
    end
    n = length(s)
    n == 1 && throw(ArgumentError("Supply both endpoints or an empty parameter vector"))
    for i in 1:n-1
        δ = (s[i+1] - s[i]) / (refinement + 1)
        for j in 0:refinement
            d.params[i, j+1] = s[i] + j * δ
        end
    end
    for i in 1:n-1, j in 1:refinement+1
        d.points[i, j] = d.path(d.params[i, j])
    end
    resize!(d.next, max(n-1, 0))
    if n > 1
        d.next[1:n-2] .= 2:n-1
        d.next[n-1] = 0
    end
    d.nref = refinement
    return d
end

"""
    add_node!(d::DiscretizedPath, idx)

Add a new node to the discretization, and return the indexes of all affected points. The indexes
are valid on the `points` and `params` fields.

# Arguments
- `d`: a DiscretizedPath object
- `idx`: a 2-element tuple, vector, or `CartesianIndex` into the `params` field. This identifies
the point to be promoted to a node.

# Returns
- A 2-element vector of `CartesianIndices` into the `params` and `points` fields.

See also [`DiscretizedPath`](@ref), [`collect`](@ref).
"""
function add_node!(d::DiscretizedPath, idx)
    if length(idx) != 2
        throw(ArgumentError("Second argument must be a 2-element tuple or CartesianIndex"))
    end
    i, j = idx[1], idx[2]
    1 <= i <= length(d.next) && 1 <= j <= d.nref+1 ||
        throw(KeyError("Indicated new node is not in the discretization"))
    j == 1 && return
    n = length(d.next) + 1
    if n == size(d.points, 1)
        throw(MaxRefinementException("Cannot add more points to this discretization"))
    end
    s_left, s_new = d.params[i, 1], d.params[i, j]
    succ = d.next[i]
    s_right = succ == 0 ? length(d.path) : d.params[succ, 1]
    δ_left = (s_new - s_left) / (1 + d.nref)
    δ_right = (s_right - s_new) / (1 + d.nref)
    if s_left + δ_left == s_left || s_new + δ_right == s_new
        throw(MaxRefinementException("Minimum refinement length reached"))
    end

    # Replace the old interval, then append its new successor.
    for j in 2:d.nref+1
        d.params[i, j] = s_left + (j - 1) * δ_left
        d.points[i, j] = point(d.path, d.params[i, j])
    end
    for j in 1:d.nref+1
        d.params[n, j] = s_new + (j - 1) * δ_right
        d.points[n, j] = point(d.path, d.params[n, j])
    end
    d.next[i] = n
    push!(d.next, succ)
    # return indexes into the new parts
    return [CartesianIndices((i:i, 1:d.nref+1)); CartesianIndices((n:n, 1:d.nref+1))]
end

"""
    collect(d::DiscretizedPath, which=:nodes)

Collect the points and parameters of a discretized path.

# Arguments
- `d`: a DiscretizedPath object
- `which`: return the nodes if :nodes, test points if :test, or all if :all

# Returns
- Tuple of two vectors: parameter values and points on the path

See also [`add_node!`](@ref), [`DiscretizedPath`](@ref).
"""
function Base.collect(d::DiscretizedPath, which=:nodes)
    n = length(d.next)
    if which == :nodes
        columns = 1:1
    elseif which == :test
        columns = 2:d.nref+1
    elseif which == :all
        columns = 1:d.nref+1
    else
        throw(ArgumentError("Second argument must be :nodes, :test, or :all"))
    end
    n == 0 && return similar(d.params, 0), similar(d.points, 0)
    params = similar(d.params, length(columns), n)
    points = similar(d.points, length(columns), n)
    k = idx = 1
    while idx != 0
        @views params[:, k] .= d.params[idx, columns]
        @views points[:, k] .= d.points[idx, columns]
        idx = d.next[idx]
        k += 1
    end
    if which in (:nodes, :all)
        s_fin = length(d.path)
        return [vec(params); s_fin], [vec(points); point(d.path, s_fin)]
    else
        return vec(params), vec(points)
    end
end

# Find the distance to the nearest neighbor at every point.
function spacing(d::DiscretizedPath{T,F}) where {T,F}
    m = length(d.next)
    n = d.nref + 1
    Δ = fill(F(Inf), m, n)
    m == 0 && return Δ
    if n == 1
        for i in eachindex(d.next)
            succ = d.next[i]
            next_point = succ == 0 ? point(d.path, length(d.path)) : d.points[succ, 1]
            δ = abs(next_point - d.points[i, 1])
            Δ[i, 1] = min(Δ[i, 1], δ)
            succ == 0 || (Δ[succ, 1] = min(Δ[succ, 1], δ))
        end
        return Δ
    end
    idx = 1
    while true
        δ = @. abs(d.points[idx, 2:n] - d.points[idx, 1:n-1])
        @. Δ[idx, 2:n-1] = min(δ[1:n-2], δ[2:n-1])
        Δ[idx, 1] = min(Δ[idx, 1], δ[1])
        idx1 = d.next[idx]
        if idx1 == 0
            Δ[idx, n] = δ[n-1]
            break
        end
        δ₊ = abs(d.points[idx1, 1] - d.points[idx, n])
        Δ[idx, n] = min(δ[n-1], δ₊)
        Δ[idx1, 1] = δ₊
        idx = idx1
    end
    return Δ
end
