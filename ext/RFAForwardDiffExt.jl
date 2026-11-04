module RFAForwardDiffExt

using RationalFunctionApproximation, ForwardDiff
const RFA = RationalFunctionApproximation

RFA.derivative(f::RFA.AbstractRationalFunction) = z -> derivative(f, z)

function RFA.derivative(f::RFA.AbstractRationalFunction, z::Number)
    # Zygote loads ForwardDiff. Keep both public methods in this extension.
    zygote = Base.get_extension(RFA, :RFAZygoteExt)
    isnothing(zygote) || return zygote.derivative_at(f, z)
    ζ = [reim(complex(z))...]
    g = ForwardDiff.gradient(x -> real(f(complex(x...))), ζ)
    return conj(complex(g...))
end

end # module
