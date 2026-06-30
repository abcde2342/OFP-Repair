import math
import mpmath
import sys
import scipy
import csv
import numpy as np
from scipy.special import lambertw
from scipy.special import iv
import matplotlib.pyplot as plt
import scipy.special
import scipy.stats
import scipy.spatial.distance
import torch
import warnings
from scipy.special._ufuncs import _iv_ratio_c as iv_ratio_c
from sympy import log as sympy_log
mpmath.mp.dps = 200
from ctypes import cdll
from scipy.special import dawsn
import alglib_ext
from scipy.special import shichi, sici
from scipy.special import fresnel, j1
from scipy.special import ndtri
from scipy.stats import norm
import ctypes
try:
    ctypes.CDLL("libopenblas.so", mode=ctypes.RTLD_GLOBAL)
except OSError:
    pass

def load_gsl():
    for path in ("/usr/local/lib/libgsl.so", "libgsl.so", "/usr/lib/x86_64-linux-gnu/libgsl.so"):
        try:
            return ctypes.CDLL(path)
        except OSError:
            continue
    raise OSError("Cannot load libgsl.so")

gsl = load_gsl()
gsl.gsl_sf_fermi_dirac_1.restype = ctypes.c_double  
gsl.gsl_sf_fermi_dirac_1.argtypes = [ctypes.c_double]  
gsl.gsl_sf_fermi_dirac_0.restype = ctypes.c_double  
gsl.gsl_sf_fermi_dirac_0.argtypes = [ctypes.c_double]
gsl.gsl_cdf_ugaussian_Pinv.restype = ctypes.c_double
gsl.gsl_cdf_ugaussian_Pinv.argtypes = [ctypes.c_double]
gsl.gsl_sf_legendre_Q1.restype = ctypes.c_double
gsl.gsl_sf_legendre_Q1.argtypes = [ctypes.c_double]
gsl.gsl_sf_bessel_j0.restype = ctypes.c_double
gsl.gsl_sf_bessel_j0.argtypes = [ctypes.c_double]
gsl.gsl_sf_Si.restype = ctypes.c_double
gsl.gsl_sf_Si.argtypes = [ctypes.c_double]
gsl.gsl_sf_bessel_y0.restype = ctypes.c_double
gsl.gsl_sf_bessel_y0.argtypes = [ctypes.c_double]
gsl.gsl_sf_psi_1.restype = ctypes.c_double
gsl.gsl_sf_psi_1.argtypes = [ctypes.c_double]
gsl.gsl_sf_fermi_dirac_m1.restype = ctypes.c_double
gsl.gsl_sf_fermi_dirac_m1.argtypes = [ctypes.c_double]
gsl.gsl_sf_erf_Q.restype = ctypes.c_double
gsl.gsl_sf_erf_Q.argtypes = [ctypes.c_double]
gsl.gsl_sf_lambert_W0.restype = ctypes.c_double
gsl.gsl_sf_lambert_W0.argtypes = [ctypes.c_double]
gsl.gsl_sf_bessel_j1.restype = ctypes.c_double 
gsl.gsl_sf_bessel_j1.argtypes = [ctypes.c_double]  
gsl.gsl_sf_hyperg_0F1.restype = ctypes.c_double
gsl.gsl_sf_hyperg_0F1.argtypes = [ctypes.c_double, ctypes.c_double]
gsl.gsl_sf_gamma_inc_Q.restype = ctypes.c_double
gsl.gsl_sf_gamma_inc_Q.argtypes = [ctypes.c_double, ctypes.c_double]

def get_tuple_dims(t):
    if isinstance(t, tuple):
        return (len(t),) + get_tuple_dims(t[0])
    return ()


def calculate_condition(a, b):
    return abs(a/(a-b)) + abs(b/(a-b))

def exp_BI(x):
    bx = scipy.special.iv(1, x)  # equal to gsl_sf_bessel_I1(x)
    con = calculate_condition(math.exp(bx), 1.0)
    return (math.exp(bx)-1.0)/x


def bJ_sin(x):
    bx = scipy.special.j0(x)
    con = calculate_condition(1.0, bx)
    return (1-bx)/math.sin(x)

def di_tan(x):
    dx = scipy.special.spence(1.0-x)  # scipy's spence is the dilogarithm (dilog) function
    val = 1/dx - 1/math.tan(x)
    con = calculate_condition(1/dx, 1/math.tan(x))
    return 1/dx-1/math.tan(x)

def log_erf(x):
    dx = scipy.special.erf(x)
    val = 1 - dx
    con1 = abs(1/math.log(1-dx))
    con2 = abs(1/math.log(1+x))
    return math.log(1-dx)/math.log(1+x)


def acos_fd(x):
    fx = gsl.gsl_sf_fermi_dirac_1(x)
    val = (math.acos(x)**2 - 3 * fx) 
    con = calculate_condition(math.acos(x)**2, 3 * fx) 
    return (math.acos(x)**2 - 3 * fx) / x


def ei(x):
    eix = gsl.gsl_cdf_ugaussian_Pinv(x/2.0 + 0.5) / math.sqrt(2.0)
    denominator = math.cos(x) - math.exp(x)
    con1 = x*((math.sqrt(math.pi) / 2.0) * math.exp(0.5 * (gsl.gsl_cdf_ugaussian_Pinv(x/2.0 + 0.5))**2))/eix
    con2 = calculate_condition(math.cos(x), math.exp(x))
    return math.sin(eix) / (math.cos(x) - math.exp(x))


def Q1_W(x):
    lx = gsl.gsl_sf_lambert_W0(x)
    q1 = gsl.gsl_sf_legendre_Q1(x)
    return (1+q1)/(lx*lx)

def bj_tan(x):
    jx = gsl.gsl_sf_bessel_j0(x)
    return (1.0 - jx)/ (x*math.tan(x))

def Si_tan(x):
    sx = gsl.gsl_sf_Si(x)
    return (sx - math.tan(x))/(x*x*x)

def by_psi(x):
    if x >= 0:
        yx = gsl.gsl_sf_bessel_y0(x)
    else:
        yx = -gsl.gsl_sf_bessel_y0(-x)
    tx = gsl.gsl_sf_psi_1(x)
    return (yx * yx - tx)


def fdm_log(x):
    fdx = gsl.gsl_sf_fermi_dirac_m1(x)
    return (2*fdx - 1) / math.log(1+x)

def eQ_sqrt(x):
    qx = gsl.gsl_sf_erf_Q(x)
    sqrt_term = math.sqrt(1 + x)
    return (2 * qx - sqrt_term) / x


def W_var(x):
    lx = gsl.gsl_sf_lambert_W0(x)
    return (lx-1) / (lx*lx - 1)

def W_log(x):
    lx = gsl.gsl_sf_lambert_W0(x)
    return (lx-1) / (lx*math.log(x) - 1)


def pow_df(x):
    val = 1.0+alglib_ext.dawsonintegral(x)
    return pow(val, 1.0/x)


def chi_ci(x):
    shix, chix = shichi(x)
    six, cix = sici(x)
    return (chix-cix)/(x*x)

def fc_bj(x):
    sx, cx = scipy.special.fresnel(x)
    val = 1.0/cx +  gsl.gsl_sf_bessel_j1(x) - math.sin(x)/(x*x)
    return (1.0/cx +  gsl.gsl_sf_bessel_j1(x) - math.sin(x)/(x*x))


def cos_x2(x):
    return (1 - math.cos(x)) / (x*x)

def exp_2(x):
    return ((math.exp(x) - 2) + math.exp(-x)) / x

def cos_sin(x):
    return (1- math.cos(x)) / math.sin(x)

def sin_sin(x):
    eps = 1e-6
    return (math.sin(x + eps) - math.sin(x))

def tan_tan(x):
    eps = 1e-6
    return (math.tan(x + eps) - math.tan(x))

def cos_cos(x):
    eps = 1e-6
    return (math.cos(x + eps) - math.cos(x))

def exp_exp(x):
    return (math.exp(x) - math.exp(-x))

def exp_1(x):
    return (math.exp(x) - 1)

def x_tan(x):
    return (1 / x - 1 / math.tan(x))

def log_log(x):
    return math.log(1 - x) / math.log(1 + x)

def log_x(x):
    return math.log((1 - x) / (1 + x))

def sqrt_exp(x):
    return math.sqrt((math.exp(2 * x) - 1) / (math.exp(x) - 1))

def sin_tan(x):
    return  (x - math.sin(x)) / (x - math.tan(x))

def exp_x(x):
    return (math.exp(x) - 1)/x

def x_x2(x):
    return ((x - 1) / ((x * x) - 1))

func_list = [
    exp_BI,
    bJ_sin,
    di_tan,
    log_erf,
    acos_fd,
    ei,
    Q1_W,  
    bj_tan,  
    Si_tan,   
    by_psi,   
    fdm_log,   
    eQ_sqrt,   
    W_var,   
    W_log,   
    pow_df,   
    chi_ci,   
    fc_bj,   
    cos_x2,
    exp_2,
    cos_sin,
    sin_sin,
    tan_tan,
    cos_cos,
    exp_exp,
    exp_1,
    x_tan,
    log_log,
    log_x,
    sqrt_exp,
    sin_tan,
    exp_x,
    x_x2,
]

func_name_list = [
    "1 exp_BI",
    "2 bJ_sin",
    "3 di_tan",
    "4 log_erf",
    "5 acos_fd",
    "6 ei",
    "7 Q1_W",
    "8 bj_tan",
    "9 Si_tan",
    "10 by_psi",
    "11 fdm_log",
    "12 eQ_sqrt",
    "13 W_var",
    "14 W_log",
    "15 pow_df",
    "16 chi_ci",
    "17 fc_bj",
    "s1 cos_x2",
    "s2 exp_2",
    "s3 cos_sin",
    "s4 sin_sin",
    "s5 tan_tan",
    "s6 cos_cos",
    "s7 exp_exp",
    "s8 exp_1",
    "s9 x_tan",
    "s10 log_log",
    "s11 log_x",
    "s12 sqrt_exp",
    "s13 sin_tan",
    "s14 exp_x",
    "s15 x_x2",
]

center_list = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    math.e,
    math.e,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1.5707963267948966192313216916398,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1.0,
]



def sci4034_1(x):
    return 1-math.exp(-0.5*x**2)

def sci4034_2(x):
    return math.sqrt(-2*math.log(1-x))

def sci3547(y):
    return ndtri((2 - y) / 2.0) / np.sqrt(2)

def sci3545_1(x):
    return 2 * (1 - norm.cdf(1 / np.sqrt(x)))


def sci3545_2(q):
    return norm.ppf(1-q / 2)

sci_extra_list = [
    ("c2 sci4034_1", 0.0, 1e-5, sci4034_1),
    ("c3 sci4034_2", 0.0, 1e-5, sci4034_2),
    ("c4 sci3547", 0.0, 1e-5, sci3547),
    ("c5 sci3545_1", 0.0, 1e-1, sci3545_1),
    ("c1 sci3545_2", 0.0, 1e-5, sci3545_2),
]

H_VALUES = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

def case_mp_value(x):
    if isinstance(x, (complex, type(mpmath.mpc(0)))):
        return mpmath.mpc(x)
    if isinstance(x, np.generic):
        return mpmath.mpf(str(float(x)))
    return mpmath.mpf(str(x))

def case_scalar_condition(fun, x, h):
    x = case_mp_value(x)
    fx = fun(x)
    mh = mpmath.mpf(str(h))
    for step in (mh, -mh):
        try:
            d = (fun(x + step) - fx) / step
            if mpmath.isfinite(abs(d)):
                value = abs(d) if abs(fx) == 0 else abs(x * d / fx)
                return float(value)
        except Exception:
            continue
    return float("nan")

def case_vector_condition(fun, xs, h):
    xs = [case_mp_value(x) for x in xs]
    fx = fun(xs)
    values = []
    mh = mpmath.mpf(str(h))
    for i, xi in enumerate(xs):
        for step in (mh, -mh):
            try:
                perturbed = list(xs)
                perturbed[i] += step
                d = (fun(perturbed) - fx) / step
                if mpmath.isfinite(abs(d)):
                    values.append(abs(d) if abs(fx) == 0 else abs(xi * d / fx))
                    break
            except Exception:
                continue
    return float(max(values)) if values else float("nan")

def case_log_softmax_component(xs, index=0):
    vals = [case_mp_value(v) for v in xs]
    m = max(vals)
    return vals[index] - m - mpmath.log(mpmath.fsum([mpmath.e ** (v - m) for v in vals]))

def case_pnorm(xs, p):
    p = mpmath.mpf(str(p))
    return mpmath.fsum([abs(case_mp_value(x)) ** p for x in xs]) ** (1 / p)

def case_softplus(x, beta):
    y = case_mp_value(beta) * case_mp_value(x)
    if y >= 0:
        return case_mp_value(x) + mpmath.log(1 + mpmath.e ** (-y)) / beta
    return mpmath.log(1 + mpmath.e ** y) / beta

def case_sigmoid(x):
    x = case_mp_value(x)
    if x >= 0:
        return 1 / (1 + mpmath.e ** (-x))
    ex = mpmath.e ** x
    return ex / (1 + ex)

def case_normalized_sinc(x):
    x = case_mp_value(x)
    if x == 0:
        return mpmath.mpf(1)
    return mpmath.sin(mpmath.pi * x) / (mpmath.pi * x)

def case_tukey_q(p, lam):
    p = mpmath.mpf(float(p))
    lam = mpmath.mpf(str(lam))
    return (p ** lam - (1 - p) ** lam) / lam

def levy_stable_issue_grid(alpha=0.93, beta=1.0):
    dist = scipy.stats.levy_stable
    old_parameterization = dist.parameterization
    try:
        dist.parameterization = "S0"
        x = np.linspace(dist.ppf(0.9, alpha, beta), dist.ppf(0.9999, alpha, beta) * 3, 1000)
        return x[x > 0]
    finally:
        dist.parameterization = old_parameterization

def levy_stable_sf_mp(x, alpha=0.93, beta=1.0, dps=100):
    from scipy import optimize
    from scipy.stats import _levy_stable as levy_impl

    with mpmath.workdps(dps):
        zeta = -beta * np.tan(np.pi * alpha / 2.0)
        x0, alpha, beta = levy_impl._nolan_round_difficult_input(
            float(x), float(alpha), float(beta), zeta, 0.005, 0.005
        )
        nolan = levy_impl.Nolan(alpha, beta, x0)
        xi = nolan.xi
        c3 = mpmath.mpf(str(nolan.c3))
        g = nolan.g
        if alpha == 1 and beta < 0 or x0 < nolan.zeta:
            return mpmath.mpf(1) - mpmath.mpf(str(levy_impl._cdf_single_value_piecewise_Z0(x0, alpha, beta)))
        if x0 == nolan.zeta:
            return mpmath.mpf("0.5") + mpmath.mpf(str(xi)) / mpmath.pi
        left_support = -xi
        right_support = np.pi / 2

        def exp_integrand(theta):
            return np.exp(-g(theta))

        with np.errstate(all="ignore"):
            if alpha > 1:
                if exp_integrand(-xi) != 0.0:
                    res = optimize.minimize(exp_integrand, (-xi,), method="L-BFGS-B", bounds=[(-xi, np.pi / 2)])
                    left_support = res.x[0]
            elif exp_integrand(np.pi / 2) != 0.0:
                res = optimize.minimize(exp_integrand, (np.pi / 2,), method="L-BFGS-B", bounds=[(-xi, np.pi / 2)])
                right_support = res.x[0]

        def sf_integrand(theta):
            return -mpmath.expm1(-mpmath.mpf(str(float(g(float(theta))))))

        return c3 * mpmath.quad(sf_integrand, [mpmath.mpf(str(left_support)), mpmath.mpf(str(right_support))])

def low_case_isfinite(value):
    try:
        return bool(np.all(np.isfinite(value)))
    except TypeError:
        return math.isfinite(abs(value))

def low_case_scalar_condition(fun, x, h):
    try:
        fx = fun(x)
    except Exception as exc:
        return float("nan"), f"error:{type(exc).__name__}:{exc}"
    if not low_case_isfinite(fx):
        return float("nan"), "nonfinite_f_x"
    for step in (h, -h):
        try:
            f_neighbor = fun(x + step)
            if not low_case_isfinite(f_neighbor):
                continue
            derivative = (f_neighbor - fx) / step
            if not low_case_isfinite(derivative):
                continue
            condition = abs(derivative) if abs(fx) == 0 else abs(x * derivative / fx)
            return float(condition), "ok"
        except Exception:
            continue
    return float("nan"), "nonfinite_or_invalid_neighbor"

def low_case_vector_condition(fun, xs, h):
    xs = list(xs)
    try:
        fx = fun(xs)
    except Exception as exc:
        return float("nan"), f"error:{type(exc).__name__}:{exc}"
    if not low_case_isfinite(fx):
        return float("nan"), "nonfinite_f_x"
    values = []
    for i, xi in enumerate(xs):
        for step in (h, -h):
            try:
                perturbed = list(xs)
                perturbed[i] = perturbed[i] + step
                f_neighbor = fun(perturbed)
                if not low_case_isfinite(f_neighbor):
                    continue
                derivative = (f_neighbor - fx) / step
                if not low_case_isfinite(derivative):
                    continue
                value = abs(derivative) if abs(fx) == 0 else abs(xi * derivative / fx)
                values.append(float(value))
                break
            except Exception:
                continue
    if not values:
        return float("nan"), "nonfinite_or_invalid_neighbor"
    return max(values), "ok"

def low_log_softmax_component(xs, index=0):
    vals = np.asarray(xs, dtype=np.float32)
    return float(scipy.special.log_softmax(vals)[index])

def scipy_log_softmax_naive(x):
    return scipy.special.log_softmax(x)

def scipy_log_softmax_component(xs, index=0):
    return float(scipy_log_softmax_naive(xs)[index])

def log_softmax_naive(x):
    x_np = np.asarray(x)
    x_torch = torch.tensor(x_np, dtype=torch.float32)
    return torch.log_softmax(x_torch, dim=0).detach().cpu().numpy()

def torch_log_softmax_component(xs, index=0):
    return float(log_softmax_naive(xs)[index])

def log_softmax(x, dim=-1):
    return x - torch.logsumexp(x, dim=dim, keepdim=True)

def log_softmax_original_component(x, index=0):
    a = torch.tensor(np.asarray(x, dtype=np.float32))
    pt = log_softmax(a, dim=-1)
    return float(pt[index].item())

def pnorm_naive(vals, p):
    arr = torch.as_tensor(vals, dtype=torch.float64)
    x1 = arr.reshape(1, -1)
    x2 = torch.zeros_like(x1)
    out = torch.nn.PairwiseDistance(p=p)(x1, x2)
    return float(out)

def low_pnorm(xs, p):
    vals = np.asarray(xs, dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        return float(np.sum(np.abs(vals) ** p) ** (1.0 / p))

def low_softplus(x, beta):
    return math.log1p(math.exp(beta * x)) / beta

def low_sigmoid(x):
    return math.exp(x) / (1.0 + math.exp(x))

def low_normalized_sinc(x):
    return float(np.sinc(x))

def pytorch_result():
    zero = torch.zeros(1, requires_grad=True)
    dx = torch.autograd.grad(
        torch.sinc(zero),
        zero,
        create_graph=True
    )[0]
    ddx = torch.autograd.grad(dx, zero)[0]
    return ddx.item()

def pytorch_sinc_second_derivative_condition(h):
    return case_scalar_condition(case_normalized_sinc, 0, h)

def low_tukey_q(p, lam):
    return scipy.stats.tukeylambda.ppf(p, lam)

def low_levy_stable_sf(x, alpha=0.93, beta=1.0):
    dist = scipy.stats.levy_stable
    old_parameterization = dist.parameterization
    try:
        dist.parameterization = "S0"
        return float(dist.sf(x, alpha, beta))
    finally:
        dist.parameterization = old_parameterization

def low_complex_from_parts(xs):
    return complex(xs[0], xs[1])

def low_hyp1f1_real_part(xs):
    z = complex(xs[2], xs[3])
    return float(scipy.special.hyp1f1(xs[0], xs[1], z).real)

def high_hyp1f1_real_part(xs):
    z = mpmath.mpc(xs[2], xs[3])
    return mpmath.hyp1f1(xs[0], xs[1], z)

def low_log1mexp(x):
    return math.log(1.0 - math.exp(x))

def low_cloglog(p):
    return math.log(-math.log(1.0 - p))

def low_log_cosh(x):
    return math.log(math.cosh(x))

def low_log_ratio_naive(b):
    return math.log((b / (b + 1.0)) / b)

def condition_iv_ratio_c(h):
    v = 0.500000000001
    x = 25.0

    def f(v, x):
        return float(iv_ratio_c(v, x))

    y = f(v, x)
    df_dv = (f(v + h, x) - y) / h
    return abs(v * df_dv / y)

def condition_kolmogi(h):
    p = 0.99999999999999

    def f_kolmogi(p):
        return float(scipy.special.kolmogi(p))

    y = f_kolmogi(p)
    df_dp = (y - f_kolmogi(p - h)) / h
    return abs(p * df_dp / y)

def condition_ksone_pdf_near_1(h):
    n = 3
    x = 1.0 - 1e-6

    def f_ksone_pdf(x):
        return float(scipy.stats.ksone(n).pdf(x))

    y = f_ksone_pdf(x)
    df_dx = (y - f_ksone_pdf(x - h)) / h
    return abs(x * df_dx / y)

def condition_ksone2_pdf_near_half(h):
    n = 2
    x = 0.5000000001

    def f_ksone2_pdf(x):
        return float(scipy.stats.ksone(n).pdf(x))

    y = f_ksone2_pdf(x)
    df_dx = (f_ksone2_pdf(x + h) - f_ksone2_pdf(x - h)) / (2 * h)
    return abs(x * df_dx / y)

def condition_smirnovi_n1(h):
    n = 1
    p = 1.0 - 1e-6

    def f_smirnovi_n1(p):
        return float(scipy.special.smirnovi(n, p))

    y = f_smirnovi_n1(p)
    df_dp = (y - f_smirnovi_n1(p - h)) / h
    return abs(p * df_dp / y)

def condition_sympy_log_close_to_1(h):
    x = 1.00000001

    def f(x):
        return sympy_log(x)

    y = f(x)
    df_dx = (f(x + h) - y) / h
    return float(abs(x * df_dx / y))

def calculate_case_conditions():
    x0 = 1 - np.log(2 * np.finfo(np.float32).eps)
    xgrid = levy_stable_issue_grid()
    scipy_17194_x = float(xgrid[333])
    cases = [
        ("case numpy/numpy#22609", lambda h: low_case_vector_condition(lambda xs: np.log1p(low_complex_from_parts(xs)), [1e-18, 1e-18], h)[0]),
        ("case numpy/numpy#5687", lambda h: low_case_vector_condition(lambda xs: 0j, [616.47292227535877, 53.814558958179042], h)[0]),
        ("case numpy/numpy#6081", lambda h: low_case_vector_condition(lambda xs: np.arccos(low_complex_from_parts(xs)).imag, [0.0, 1e-14], h)[0]),
        ("case numpy/numpy#6082", lambda h: low_case_vector_condition(lambda xs: np.arcsin(low_complex_from_parts(xs)).imag, [0.01, 1e-14], h)[0]),
        ("case numpy/numpy#6083", lambda h: low_case_vector_condition(lambda xs: np.arctan(low_complex_from_parts(xs)).imag, [0.01, 1e-14], h)[0]),
        ("case pytorch/pytorch#113708", lambda h: low_case_vector_condition(lambda xs: torch_log_softmax_component(xs, 0), [x0, 0.0], h)[0]),
        ("case pytorch/pytorch#171249", lambda h: case_scalar_condition(lambda x: case_softplus(x, 1e30), 2.0, h)),
        ("case pytorch/pytorch#184036", lambda h: low_case_vector_condition(lambda xs: pnorm_naive(xs, 100), [2.5577, 2.5867, 2.5688, 3.2271], h)[0]),
        ("case pytorch/pytorch#39242", lambda h: low_case_scalar_condition(low_log1mexp, -1e-16, h)[0]),
        ("case pytorch/pytorch#56340", lambda h: low_case_vector_condition(lambda xs: log_softmax_original_component(xs, 0), [800.0, 800.0, 400.0], h)[0]),
        ("case pytorch/pytorch#89459", pytorch_sinc_second_derivative_condition),
        ("case scipy/scipy#17194", lambda h: low_case_scalar_condition(low_levy_stable_sf, scipy_17194_x, h)[0]),
        ("case scipy/scipy#19521", lambda h: low_case_vector_condition(lambda xs: scipy_log_softmax_component(xs, 0), [x0, 0.0], h)[0]),
        ("case scipy/scipy#21114", lambda h: case_vector_condition(lambda xs: case_pnorm(xs, 1_000_000), [10.0, 20.0, 3.3, 7.3, 8.8], h)),
        ("case scipy/scipy#21370", lambda h: low_case_scalar_condition(lambda p: low_tukey_q(p, 0.2), 0.501, h)[0]),
        ("case scipy/scipy#5349", lambda h: case_vector_condition(high_hyp1f1_real_part, [30, 70, 0, 20], h)),
        ("case statsmodels/statsmodels#1604", lambda h: case_scalar_condition(case_sigmoid, mpmath.mpf("1741.90624"), h)),
        ("case statsmodels/statsmodels#826", lambda h: low_case_scalar_condition(low_cloglog, 1e-10, h)[0]),
        ("case sympy/sympy#12671", lambda h: case_scalar_condition(lambda x: mpmath.log(mpmath.cosh(x)), 1000.0, h)),
        ("case sympy/sympy#28514", lambda h: low_case_scalar_condition(low_log_ratio_naive, 1e-18, h)[0]),
        ("case 1: iv_ratio_c", condition_iv_ratio_c),
        ("case 2: kolmogi", condition_kolmogi),
        ("case 3: ksone.pdf near 1", condition_ksone_pdf_near_1),
        ("case 4: ksone(2).pdf near 0.5", condition_ksone2_pdf_near_half),
        ("case 5: smirnovi(1, p), p close to 1", condition_smirnovi_n1),
        ("case 6: log(x), x close to 1", condition_sympy_log_close_to_1),
    ]
    rows = []
    for name, compute in cases:
        values = []
        for h in H_VALUES:
            try:
                values.append(compute(h))
            except Exception as exc:
                values.append(float("nan"))
        rows.append((name, "", *values))
    return rows

def calculate_con(x, fun, h):
    y = fun(x)
    der = (fun(x+h)-y) / h
    con = x * der / y
    return round(abs(con), 10)

def calculate_func_list_conditions():
    rows = []
    x_offset = 1e-2
    for name, center, fun in zip(func_name_list, center_list, func_list):
        x = center + x_offset
        values = []
        for h in H_VALUES:
            try:
                values.append(calculate_con(x, fun, h))
            except Exception as exc:
                values.append(float("nan"))
        rows.append((name, x, *values))
    for name, center, x, fun in sci_extra_list:
        values = []
        for h in H_VALUES:
            try:
                values.append(calculate_con(x, fun, h))
            except Exception as exc:
                values.append(float("nan"))
        rows.append((name, x, *values))
    rows.extend(calculate_case_conditions())
    return rows

def h_column_name(h):
    return f"h_{h:.0e}".replace("+", "")

def csv_number(value):
    if value == "":
        return value
    try:
        value = float(value)
    except (TypeError, ValueError):
        return value
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return f"{value:.10e}"

def write_func_list_conditions(path="condition_number_func_list.csv"):
    rows = calculate_func_list_conditions()
    with open(path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["function_id", "x", *[h_column_name(h) for h in H_VALUES]])
        for row in rows:
            writer.writerow([row[0], *[csv_number(value) for value in row[1:]]])
    return rows

if __name__ == "__main__":
    rows = write_func_list_conditions()
    for row in rows:
        print(row[0], row[2:])
