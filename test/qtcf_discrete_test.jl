# Run the shared experiments on fixed discrete grids.
const QTCF_EXPERIMENT_MODE = :discrete
include("qtcf_test.jl")
abspath(PROGRAM_FILE) == (@__FILE__) && run_experiments()
