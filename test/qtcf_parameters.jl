# Fitting parameters.
const T = Float64
const METHODS = (:AAA, :TCF, :QTCF)
const EXPERIMENT_EVAL_METHOD = OneDiv()
const EXPERIMENT_WEIGHT_METHODS = (TCF=OneDiv(), QTCF=Classic())
const EXPERIMENT_PARAMETERS = (
    interval=(float_type=T, tol=T(1e-11), allowed=true, max_degree=200,
        initial_refinement=15, refinement=4, stagnation=50),
    schwarz=(float_type=T, tol=T(4e-11), allowed=true, max_degree=200,
        initial_refinement=28, refinement=9, stagnation=50),
)
const AAA_LOEWNER_TOL = 100eps(T)
# Fixed-point startup on the normalized fitting curve.
const QTCF_START_PARAMETER = 0.5

# All paper curves are paired by reflection of the normalized parameter.
qtcf_parameter_reflection(s) = 1.0 - s

function experiment_parameters(family, method; discrete=false)
    p = EXPERIMENT_PARAMETERS[family]
    max_iter = (method === :TCF ? 2 : 1) * p.max_degree
    options = (; p.float_type, p.tol, p.allowed, p.stagnation, max_iter)
    if discrete
        method === :AAA && return (; options..., loewner_tol=AAA_LOEWNER_TOL)
        return method === :QTCF ? (; options..., p.max_degree) : options
    end
    options = (; options..., p.max_degree, p.refinement)
    return method === :QTCF ? (; options..., initial_refinement=p.initial_refinement) : options
end

# Sampling and validation.
const DISCRETE_UNIFORM_N = 5000
const DISCRETE_LOG_EXPONENTS = -7:-3
const ABS_UNIFORM_N = 10001
const ABS_LOG_EXPONENTS = -0.1 .* collect(10:1000)
const ABS_FITTING_LOG_EXPONENTS = -52:1//10:-1
const INTERVAL_DISPLAY_DEGREE = 100
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
const POLE_RESIDUE_LOG_RANGE = (-4.0, 0.0)
const FIGURE_TITLES = (
    ("Figure_1", "Ellipse (1, 1/2)"),
    ("Figure_6_1_interval_convergence", "Interval convergence"),
    ("Figure_6_2_continuum_poles_1", "First four Schwarz curves"),
    ("Figure_6_3_continuum_poles_2", "Remaining Schwarz curves"),
)
