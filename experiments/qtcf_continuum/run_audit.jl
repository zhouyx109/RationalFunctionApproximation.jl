# Read the completed fits and run the fixed-node normalization checks.
include("certify_intervals.jl")
include("audit_update.jl")
audit_update(deserialize(joinpath(DEST, "results.jls")))
