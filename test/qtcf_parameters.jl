# Fitting parameters.
const T = Float64
const METHODS = (:AAA, :TCF, :QTCF)
const EXPERIMENT_PARAMETERS = (
    interval=(float_type=T, tol=100eps(T), allowed=true, max_degree=200,
        initial_refinement=28, refinement=4, stagnation=200),
    schwarz=(float_type=T, tol=1e-12, allowed=true, max_degree=200,
        initial_refinement=15, refinement=4, stagnation=50),
)
const AAA_LOEWNER_TOL = 100eps(T)

function experiment_parameters(family, method; discrete=false)
    p = EXPERIMENT_PARAMETERS[family]
    options = (; (key => value for (key, value) in pairs(p)
        if key !== :initial_refinement && (!discrete || key ∉ (:max_degree, :refinement)))...,
        max_iter=(method === :TCF ? 2 : 1) * p.max_degree)
    discrete && return method === :AAA ? (; options..., loewner_tol=AAA_LOEWNER_TOL) : options
    return method === :QTCF ? (; options..., initial_refinement=p.initial_refinement) : options
end

# Sampling and validation.
const DISCRETE_UNIFORM_N = 5000
const DISCRETE_LOG_MIN, DISCRETE_LOG_STEP, DISCRETE_LOG_MAX = -16, 0.1, -3
const DISCRETE_LOG_EXPONENTS = DISCRETE_LOG_MIN:DISCRETE_LOG_STEP:DISCRETE_LOG_MAX
const ABS_UNIFORM_N = 10001
const ABS_LOG_EXPONENTS = -100:0.1:-1
const SCHWARZ_AUDIT_N = 1800
const SCHWARZ_CURVE_N = 700
const INLET_OFFSET = 0.1
const INLET_COS = (-0.1, -0.7, -0.1)
const INLET_SIN = (-0.8, -0.4, -0.2)

# Figure labels and display.
const AAA_COLOR = RGBf(0 / 255, 114 / 255, 178 / 255)
const TCF_COLOR = RGBf(230 / 255, 159 / 255, 0 / 255)
const QTCF_COLOR = RGBf(0 / 255, 158 / 255, 115 / 255)
const CASE_LABELS = Dict(
    "ellipse_1_025" => "Ellipse (1, 1/4)", "ellipse_1_05" => "Ellipse (1, 1/2)",
    "upper_ellipse_1_05" => "Upper ellipse (1, 1/2)", "five_lobed" => "Five-lobed curve",
    "smooth_s" => "Smooth S", "two_semicircle_s" => "Two-semicircle S",
    "superellipse_6" => "Superellipse", "analytic_inlet" => "Analytic inlet",
)
const CONTINUUM_LIMITS = (-2.25, 2.25, -2.25, 2.25)
const CONTINUUM_TICKS = -2:1:2
const FIGURE_TITLES = (
    ("Figure_1", "Ellipse (1, 1/2)"),
    ("Figure_6_1_interval_convergence", "Interval convergence"),
    ("Figure_6_2_continuum_poles_1", "First four Schwarz curves"),
    ("Figure_6_3_continuum_poles_2", "Remaining Schwarz curves"),
)
