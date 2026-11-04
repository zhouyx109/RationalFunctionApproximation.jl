module RFAZygoteExt

using RationalFunctionApproximation, Zygote
const RFA = RationalFunctionApproximation

function derivative_at(f::RFA.AbstractRationalFunction, z::Number)
    ζ = complex(z)
    slope = gradient(real ∘ f, ζ)[1]
    return isnothing(slope) ? zero(ζ) : conj(slope)
end

end # module
