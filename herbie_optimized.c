#include <math.h>

double herbie_cos_x2(double x) {
    return fma(x * x, fma(x * x, 0.001388888888888889, -0.041666666666666664), 0.5);
}

double herbie_exp_2(double x) {
    return fma(x * x, x * fma(x, x * 0.002777777777777778, 0.08333333333333333), x);
}

double herbie_cos_sin(double x) {
    return x * fma(x * x, fma(x, x * 0.004166666666666667, 0.041666666666666664), 0.5);
}

double herbie_sin_sin(double x) {
    return (cos(5e-7) * cos(x) - sin(5e-7) * sin(x)) * (sin(5e-7) * 2.0);
}

double herbie_tan_tan(double x) {
    return fma(sin(x), fma(-0.5, x * x, -1.0), tan(fma(x, x, -1e-12) / (x - 1e-6)));
}

double herbie_cos_cos(double x) {
    return sin(5e-7) * sin(fma(fma(x, x, -1e-12), 1.0 / (x - 1e-6), x) * 0.5) * -2.0;
}

double herbie_exp_exp(double x) {
    return fma(x * x * x, 0.3333333333333333, x * 2.0);
}

double herbie_exp_1(double x) {
    return fma(fma(x, fma(x, 0.041666666666666664, 0.16666666666666666), 0.5), x * x, x);
}

double herbie_x_tan(double x) {
    return x / fma(x, x * -0.2, 3.0);
}

double herbie_log_log(double x) {
    return x * fma(x, fma(x, fma(x, -0.25, -0.3333333333333333), -0.5), -1.0) /
           fma(fma(x, fma(x, -0.25, 0.3333333333333333), -0.5), x * x, x);
}

double herbie_log_x(double x) {
    return fma(x * x * x, -0.6666666666666666, x * -2.0);
}

double herbie_sqrt_exp(double x) {
    return sqrt(exp(x) + 1.0);
}

double herbie_sin_tan(double x) {
    return fma(0.225, x * x, -0.5);
}

double herbie_exp_x(double x) {
    return fma(x, fma(x, fma(x, 0.041666666666666664, 0.16666666666666666), 0.5), 1.0);
}

double herbie_x_x2(double x) {
    return (x - 1.0) / fma(x, x, -1.0);
}

double herbie_sci4034_1(double x) {
    return x * fma(x * x * x, -0.125, x * 0.5);
}

double herbie_sci4034_2(double x) {
    return sqrt(-2.0 * log1p(-x));
}

double herbie_sympy_28514_log_ratio(double b) {
    return b * fma(b, fma(b, fma(b, 0.25, -0.3333333333333333), 0.5), -1.0);
}

double herbie_sympy_12671_log_cosh(double x) {
    return log(0.001388888888888889) + fma(6.0, log(x), 30.0 / (x * x));
}

double herbie_statsmodels_1604_logistic(double eta) {
    (void)eta;
    return 0.5;
}

double herbie_statsmodels_826_cloglog(double p) {
    return fma(p, fma(p, fma(p, 0.125, 0.20833333333333334), 0.5), log(p));
}

double herbie_scipy_17194_expm1_kernel(double g) {
    return fma(fma(g, fma(g, -0.041666666666666664, 0.16666666666666666), -0.5), g * g, g);
}

double herbie_pytorch_39242_log1mexp(double x) {
    return log(-expm1(x));
}
