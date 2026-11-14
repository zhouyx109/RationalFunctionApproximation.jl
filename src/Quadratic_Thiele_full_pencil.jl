# Assemble H_m for numerator roots or its trailing submatrix K_m for poles.
# Both have one row per block. Each link satisfies -ℓ_j(t)υ_j(t) = a_j(t).
function _qtcf_block_pencil(r::QuadraticThiele{T,S}, first_block::Integer) where {T,S}
    pairs, bpolys, apolys = r.pairs, r.bpolys, r.apolys
    n = max(length(bpolys) - first_block + 1, 0)
    n == 0 && return zeros(S, 0, 0), zeros(S, 0, 0)
    c_diag, d_diag = zeros(S, n), zeros(S, n)
    c_lower, d_lower = zeros(S, n - 1), ones(S, n - 1)
    c_upper, d_upper = zeros(S, n - 1), zeros(S, n - 1)
    for (j, k) in enumerate(first_block:length(bpolys))
        b = bpolys[k]
        c_diag[j] = b[1]
        length(b) == 2 && (d_diag[j] = b[2])
        j == n && continue

        z₁, z₂ = pairs[k]
        a = apolys[k]
        multiplier = last(a)
        c_lower[j] = -z₁
        if length(a) == 2
            c_upper[j] = -multiplier
        else
            c_upper[j], d_upper[j] = multiplier * z₂, -multiplier
        end
    end
    C = diagm(0 => c_diag, -1 => c_lower, 1 => c_upper)
    D = diagm(0 => d_diag, -1 => d_lower, 1 => d_upper)
    return C, D
end
