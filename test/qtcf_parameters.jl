# Shared presets, validation grids, and curve definitions.
const T = Float64
const METHODS = (:AAA, :TCF, :QTCF)
const AAA_COLOR = RGBf(0 / 255, 114 / 255, 178 / 255)
const TCF_COLOR = RGBf(230 / 255, 159 / 255, 0 / 255)
const QTCF_COLOR = RGBf(0 / 255, 158 / 255, 115 / 255)
const DISCRETE_UNIFORM_N = 5000
const DISCRETE_LOG_MIN = -16
const DISCRETE_LOG_STEP = 0.1
const DISCRETE_LOG_MAX = -3
const DISCRETE_LOG_EXPONENTS = DISCRETE_LOG_MIN:DISCRETE_LOG_STEP:DISCRETE_LOG_MAX
const SCHWARZ_AUDIT_N = 1800
const SCHWARZ_CURVE_N = 700
const RFA_INITIAL_REFINEMENT = 15

# The package converts the degree limit to each method's node/block limit.
const CONTINUUM_PARAMETERS = (
    float_type=T, tol=1e-11, allowed=true, max_degree=200,
    refinement=4, stagnation=50,
)
const CONTINUUM_QTCF_PARAMETERS = (
    ; CONTINUUM_PARAMETERS..., initial_refinement=RFA_INITIAL_REFINEMENT)

const INTERVAL_TCF_PARAMETERS = EXPERIMENT_MODE === :continuum ? CONTINUUM_PARAMETERS : (
    tol=100eps(T), allowed=true, max_iter=241, refinement=3, stagnation=241,
)
const INTERVAL_QTCF_PARAMETERS = EXPERIMENT_MODE === :continuum ? CONTINUUM_QTCF_PARAMETERS : (
    tol=100eps(T), allowed=true, max_iter=150,
    refinement=4, initial_refinement=28, stagnation=170,
)

const SCHWARZ_PARAMETERS = EXPERIMENT_MODE === :continuum ? (
    AAA=CONTINUUM_PARAMETERS, TCF=CONTINUUM_PARAMETERS,
    QTCF=CONTINUUM_QTCF_PARAMETERS,
) : let common = (tol=1e-12, allowed=true, stagnation=10)
    (AAA=(; common..., max_iter=150, refinement=4),
        TCF=(; common..., max_iter=500, refinement=3),
        QTCF=(; common..., max_iter=100, refinement=3, initial_refinement=20))
end

# Independent interval validation uses uniform and logarithmic samples.
fine = T(2) .^ (-100:0.1:-1)
abs_check = sort!(unique(vcat(range(-one(T), one(T), 10001),
    fine, -fine, fine .- one(T))))

interval_cases = (
    (
        slug="real",
        name="|x| on [-1,1]",
        title=raw"$|x|$ on $[-1,1]$",
        f=x -> abs(real(x)),
        domain=Segment{T}(-1, 1),
        points=abs_check,
    ),
    (
        slug="imaginary",
        name="|z| on [-i,i]",
        title=raw"$|z|$ on $[-i,i]$",
        f=z -> abs(z),
        domain=Segment(-1im, 1im),
        points=im .* abs_check,
    ),
)

# Pi-based curves retain exact fixed points of the supplied symmetries.
half_turn(c) = z -> 2c - z
line_reflection(c, factor) = z -> begin
    w = z - c
    c + factor * complex(real(w), -imag(w))
end

upper_ellipse_point(t) = begin
    cospi(1 - t) + 0.5im * sinpi(1 - t)
end
upper_ellipse_deriv(t) = begin
    pi * sinpi(1 - t) - 0.5im * pi * cospi(1 - t)
end
paper_upper_ellipse = Curve(upper_ellipse_point, upper_ellipse_deriv, 0, 1)

paper_lobed_full = ClosedCurve(
    t -> (1 + 0.2sinpi(5 * (t / pi))) * cispi(t / pi),
    t -> cispi(t / pi) * (cospi(5 * (t / pi)) + im * (1 + 0.2sinpi(5 * (t / pi)))),
    0, 2pi)

paper_smooth_s_point(t) = -0.75sinpi(2t) + im * (2 - 4t)
paper_smooth_s_deriv(t) = -1.5pi * cospi(2t) - 4im
paper_smooth_s_curve = Curve(paper_smooth_s_point, paper_smooth_s_deriv, 0, 1)

function paper_semicircle_s_point(t)
    if t <= 0.5
        s = 2t
        return im + cispi(0.5 + s)
    end
    s = 2t - 1
    -im + cispi(0.5 - s)
end

paper_semicircle_s_deriv(t) =
    t <= 0.5 ? 2pi * im * cispi(0.5 + 2t) : -2pi * im * cispi(0.5 - (2t - 1))
paper_semicircle_s_curve =
    Curve(paper_semicircle_s_point, paper_semicircle_s_deriv, 0, 1)

signedpow(u, p) = abs(u) < 1e-14 ? 0.0 : sign(u) * abs(u)^p
super6_point(t) = if EXPERIMENT_MODE === :discrete
    signedpow(cos(t), 1 / 3) + im * signedpow(sin(t), 1 / 3)
else
    c, s = cospi(2t), sinpi(2t)
    cispi(2t) / (c^6 + s^6)^(1 / 6)
end
super6_deriv(t) = if EXPERIMENT_MODE === :discrete
    (super6_point(t + 1e-6) - super6_point(t - 1e-6)) / 2e-6
else
    c, s = cospi(2t), sinpi(2t)
    g = c^6 + s^6
    2pi * cispi(2t) * (g^(-7 / 6) * (s * c^5 - c * s^5) + im * g^(-1 / 6))
end
paper_superellipse6 = ClosedCurve(super6_point, super6_deriv, 0,
    EXPERIMENT_MODE === :continuum ? 1 : 2pi)

const inlet_a0 = 0.1
const inlet_A = [-0.1, -0.7, -0.1]
const inlet_B = [-0.8, -0.4, -0.2]

paper_inlet_point(t) = inlet_a0 + sum(
    inlet_A[k] * cospi(k * t / pi) + im * inlet_B[k] * sinpi(k * t / pi)
    for k in eachindex(inlet_A))
paper_inlet_deriv(t) = sum(
    -k * inlet_A[k] * sinpi(k * t / pi) + im * k * inlet_B[k] * cospi(k * t / pi)
    for k in eachindex(inlet_A))
paper_inlet = ClosedCurve(paper_inlet_point, paper_inlet_deriv, 0, 2pi)

origin = 0.0 + 0.0im

qtcf_cases = NamedTuple{(:slug, :name, :domain, :reflection,
    :function_reflection, :formula, :rotation)}.([
    ("ellipse_1_025", "ellipse (1, 1/4) Schwarz", Shapes.ellipse(1, 0.25),
        conj, conj, raw"$z(t)=\cos t+\frac{i}{4}\sin t$", 1.0 + 0.0im),
    ("ellipse_1_05", "ellipse (1, 1/2) Schwarz", Shapes.ellipse(1, 0.5),
        conj, conj, raw"$z(t)=\cos t+\frac{i}{2}\sin t$", 1.0 + 0.0im),
    ("upper_ellipse_1_05", "upper ellipse (1, 1/2) Schwarz",
        paper_upper_ellipse, line_reflection(origin, -1.0 + 0.0im),
        y -> -conj(y),
        raw"$z(t)=\cos(\pi(1-t))+\frac{i}{2}\sin(\pi(1-t))$", 1im),
    ("five_lobed", "five-lobed Schwarz", paper_lobed_full,
        line_reflection(origin, -1.0 + 0.0im), y -> -conj(y),
        raw"$z(t)=(1+0.2\sin 5t)e^{it}$", 1im),
    ("smooth_s", "smooth letter S Schwarz", paper_smooth_s_curve,
        half_turn(origin), y -> -y, raw"$z(t)=-0.75\sin(2\pi t)+i(2-4t)$", 1.0 + 0.0im),
    ("two_semicircle_s", "two-semicircle S Schwarz", paper_semicircle_s_curve,
        half_turn(origin), y -> -y,
        raw"$z(t)=i+e^{i(\pi/2+2\pi t)}\ (0\leq t\leq1/2),\quad z(t)=-i+e^{i(3\pi/2-2\pi t)}\ (1/2<t\leq1)$", 1.0 + 0.0im),
    ("superellipse_6", "superellipse x^6+y^6=1 Schwarz", paper_superellipse6,
        conj, conj, EXPERIMENT_MODE === :continuum ?
            raw"$z(t)=e^{2\pi it}/(\cos^6(2\pi t)+\sin^6(2\pi t))^{1/6}$" :
            raw"$z(t)=\mathrm{sgn}(\cos t)|\cos t|^{1/3}+i\mathrm{sgn}(\sin t)|\sin t|^{1/3}$", 1.0 + 0.0im),
    ("analytic_inlet", "analytic right inlet Schwarz", paper_inlet,
        conj, conj,
        raw"$z(t)=0.1-0.1\cos t-0.7\cos2t-0.1\cos3t-i(0.8\sin t+0.4\sin2t+0.2\sin3t)$", 1.0 + 0.0im),
])

const CASE_LABELS = Dict(
    "ellipse_1_05" => "Ellipse (1, 1/2)",
    "ellipse_1_025" => "Ellipse (1, 1/4)",
    "upper_ellipse_1_05" => "Upper ellipse (1, 1/2)",
    "five_lobed" => "Five-lobed curve",
    "smooth_s" => "Smooth S",
    "two_semicircle_s" => "Two-semicircle S",
    "superellipse_6" => "Superellipse",
    "analytic_inlet" => "Analytic inlet",
)

const CONTINUUM_LIMITS = (-2.25, 2.25, -2.25, 2.25)
const CONTINUUM_TICKS = -2:1:2

const FIGURE_TITLES = (
    ("Figure_1", "Ellipse (1, 1/2)"),
    ("Figure_6_1_interval_convergence", "Interval convergence"),
    ("Figure_6_2_continuum_poles_1", "First four Schwarz curves"),
    ("Figure_6_3_continuum_poles_2", "Remaining Schwarz curves"),
)

discrete_parameters(parameters) = (; (key => value for (key, value) in pairs(parameters)
    if key ∉ (:refinement, :initial_refinement))...)
