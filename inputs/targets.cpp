#include <math.h>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include "target.h"
#include <libalglib/specialfunctions.h>
#include <boost/math/distributions/normal.hpp>

namespace {
int gExpansionTerms = 10;

int expansion_term_count(int firstDegree, int degreeStep, int maxTerms) {
    if (gExpansionTerms < firstDegree) {
        return 0;
    }
    int terms = (gExpansionTerms - firstDegree) / degreeStep + 1;
    return terms < maxTerms ? terms : maxTerms;
}

int ofp_series_terms(int maxTerms, int availableTerms) {
    if (maxTerms < 1) {
        return 0;
    }
    return maxTerms < availableTerms ? maxTerms : availableTerms;
}

double horner_series(double x, const double* coef, int n) {
    if (n <= 0) {
        return 0.0;
    }
    double y = coef[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        y = y * x + coef[i];
    }
    return y;
}

double even_series(double x, const double* coef, int n) {
    return horner_series(x * x, coef, n);
}

double ordinary_series(double x, const double* coef, int maxTerms) {
    return horner_series(x, coef, expansion_term_count(0, 1, maxTerms));
}

double even_degree_series(double x, const double* coef, int maxTerms) {
    return even_series(x, coef, expansion_term_count(0, 2, maxTerms));
}

double odd_degree_series(double x, const double* coef, int maxTerms) {
    return x * even_series(x, coef, expansion_term_count(1, 2, maxTerms));
}

double shifted_series(double x, double center, const double* coef, int maxTerms) {
    return horner_series(x - center, coef, expansion_term_count(0, 1, maxTerms));
}

double log1p_expansion(double x, int maxTerms) {
    double total = 0.0;
    double power = 1.0;
    int terms = expansion_term_count(1, 1, maxTerms);
    for (int n = 1; n <= terms; ++n) {
        power *= x;
        total += (n % 2 ? 1.0 : -1.0) * power / n;
    }
    return total;
}

double expm1_expansion(double x, int maxTerms) {
    double total = 0.0;
    double term = 1.0;
    int terms = expansion_term_count(1, 1, maxTerms);
    for (int n = 1; n <= terms; ++n) {
        term *= x / n;
        total += term;
    }
    return total;
}

std::complex<double> complex_log1p_expansion(std::complex<double> z, int maxTerms) {
    std::complex<double> total = 0.0;
    std::complex<double> power = 1.0;
    int terms = expansion_term_count(1, 1, maxTerms);
    for (int n = 1; n <= terms; ++n) {
        power *= z;
        total += (n % 2 ? 1.0 : -1.0) * power / (double)n;
    }
    return total;
}

double binomial_one_plus(double s, double alpha, int maxTerms) {
    double total = 1.0;
    double coeff = 1.0;
    double power = 1.0;
    int terms = expansion_term_count(1, 1, maxTerms);
    for (int n = 1; n <= terms; ++n) {
        coeff *= (alpha - (n - 1)) / n;
        power *= s;
        total += coeff * power;
    }
    return total;
}

double stable_pnorm(const std::vector<double>& vals, double p) {
    double scale = 0.0;
    for (double v : vals) {
        scale = std::max(scale, std::fabs(v));
    }
    if (scale == 0.0) {
        return 0.0;
    }
    double s = 0.0;
    bool skipped = false;
    for (double v : vals) {
        double r = std::fabs(v) / scale;
        if (!skipped && r == 1.0) {
            skipped = true;
            continue;
        }
        s += std::pow(r, p);
    }
    return scale * binomial_one_plus(s, 1.0 / p, 12);
}

std::complex<double> hyp1f1_expansion_complex(double a, double b, std::complex<double> z, int maxTerms) {
    std::complex<double> term = 1.0;
    std::complex<double> total = term;
    int terms = ofp_series_terms(gExpansionTerms, maxTerms);
    for (int n = 1; n <= terms; ++n) {
        term *= (a + n - 1.0) / (b + n - 1.0) * z / (double)n;
        total += term;
    }
    return total;
}

double tukey_q_expansion(double p, double lam, int maxTerms) {
    double u = 2.0 * (p - 0.5);
    double coeff = std::pow(2.0, 1.0 - lam);
    double sum = 0.0;
    int terms = ofp_series_terms(gExpansionTerms, maxTerms);
    for (int k = 0; k < terms; ++k) {
        int n = 2 * k + 1;
        double prod = 1.0;
        for (int j = 0; j < n; ++j) {
            prod *= lam - j;
        }
        sum += prod * std::pow(u, n) / std::tgamma(n + 1.0);
    }
    return coeff * sum / lam;
}

double polygamma_even_at_one(int order) {
    switch (order) {
        case 2:  return -2.4041138063191885708;       // -2! * zeta(3)
        case 4:  return -24.886266123440878231;       // -4! * zeta(5)
        case 6:  return -726.01147971498443532;       // -6! * zeta(7)
        case 8:  return -40419.524592934859294;       // -8! * zeta(9)
        case 10: return -3628799.5059087511546;       // -10! * zeta(11)
        default: return 0.0;
    }
}

namespace detail {
constexpr double kPi =
    3.141592653589793238462643383279502884;

constexpr double kInvSqrtPi =
    0.564189583547756286948079451560772586;

double erfc_inv_small_y(double y, int terms) {
    if (y <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (y >= 2.0) {
        return -std::numeric_limits<double>::infinity();
    }

    const double L1 = std::log(2.0 / kPi) - 2.0 * std::log(y);
    const double L2 = std::log(L1);

    const double r  = 1.0 / L1;
    const double r2 = r * r;
    const double r3 = r2 * r;

    double W = L1 - L2;
    if (terms >= 2) {
        W += L2 * r;
    }
    if (terms >= 3) {
        W += L2 * (-2.0 + L2) * r2 / 2.0;
    }
    if (terms >= 4) {
        W += L2 * (6.0 - 9.0 * L2 + 2.0 * L2 * L2) * r3 / 6.0;
    }

    double x = std::sqrt(0.5 * W);

    for (int i = 0; i < 3; ++i) {
        const double f  = std::erfc(x) - y;
        const double fp = -2.0 * kInvSqrtPi * std::exp(-x * x);

        if (fp == 0.0 || !std::isfinite(fp)) {
            break;
        }

        const double dx = f / fp;
        x -= dx;

        if (std::abs(dx) <= 4.0 * std::numeric_limits<double>::epsilon() * std::abs(x)) {
            break;
        }
    }

    return x;
}
}
}

extern "C" void set_expansion_terms(int terms) {
    gExpansionTerms = terms;
}

extern "C" int get_expansion_terms() {
    return gExpansionTerms;
}

double cos_x2(double x) {
    double val = (1 - cos(x)) / (x*x);
    return val;
}

double exp_2(double x) {
    double val = ((exp(x) - 2) + exp(-x))/x;
    return val;
}

double cos_sin(double x) {
    double val = (1 - cos(x)) / sin(x);
    return val;
}

double sin_sin(double x) {
    double eps = 1e-6;
    double val = sin(x + eps) - sin(x);
    return val;
}

double tan_tan(double x) {
    double eps = 1e-6;
    double val = tan(x + eps) - tan(x);
    return val;
}

double cos_cos(double x) {
    double eps = 1e-6;
    double val = cos(x + eps) - cos(x);
    return val;
}

double exp_exp(double x) {
    double val = exp(x) - exp(-x);
    return val;
}

double exp_1(double x) {
    double val = exp(x) - 1;
    return val;
}

double x_tan(double x) {
    double val = 1 / x - 1 / tan(x);
    return val;
}

double log_log(double x) {
    double val = log(1 - x) / log(1 + x);
    return val;
}

double log_x(double x) {
    double val = log((1 - x) / (1 + x));
    return val;
}

double sqrt_exp(double x) {
    double val = sqrt((exp(2 * x) - 1) / (exp(x) - 1));
    return val;
}

double sin_tan(double x) {
    double val = (x - sin(x)) / (x - tan(x));
    return val;
}

double exp_x(double x) {
    double val = ((exp(x) - 1) / x);
    return val;
}

double x_x2(double x) {
    double val = ((x - 1) / ((x * x) - 1));
    return val;
}

double exp_BI(double x) {
    double bx = gsl_sf_bessel_I1(x);
    double val = ((exp(bx) - 1.0) / x);
    return val;
}

double bJ_sin(double x) {
    double bx = gsl_sf_bessel_J0(x);
    double val = (1 - bx) / sin(x);
    return val;
}

double di_tan(double x) {
    double dx = gsl_sf_dilog(x);
    double val = 1/dx - 1/tan(x);
    return val;
}

double log_erf(double x) {
    double ex = gsl_sf_erf(x);
    double val = log(1-ex) / log(1+x);
    return val;
    // return (x * (-(2.0 / std::sqrt(M_PI)) + x * (-(2.0 / M_PI) + x * ((2.0 * (M_PI - 4.0)) / (3.0 * std::pow(M_PI, 1.5)) + x * ((4.0 * (M_PI - 3.0)) / (3.0 * M_PI * M_PI) + x * ((-96.0 + 40.0 * M_PI - 3.0 * M_PI * M_PI) / (15.0 * pow(M_PI, 2.5)) + x * (-4.0 * (120.0 - 60.0 * M_PI + 7.0 * M_PI * M_PI) / (45.0 * M_PI * M_PI * M_PI)))))))) / (x * (1 + x * (-1.0/2 + x * (1.0/3 + x * (-1.0/4 + x * (1.0/5 + x * (-1.0/6)))))));
}

double acos_fd(double x) {
    double fx = gsl_sf_fermi_dirac_1(x);
    double val = (acos(x)*acos(x) - 3*fx) / x;
    return val;
}

double sci3545_2(double q) {
    double val = gsl_cdf_ugaussian_Pinv(1.0 - q / 2.0);
    return val;
}

double sci4034_1(double x) {
    return 1-exp(-0.5*x*x);
}

double sci4034_2(double x) {
    return sqrt(-2*log(1-x));
}

double sci3547(double y) {
    if (y <= 0.0 || y >= 2.0) {
        return std::numeric_limits<double>::infinity();
    }
    boost::math::normal_distribution<double> normal;
    try {
        return boost::math::quantile(normal, (2.0 - y) / 2.0) / std::sqrt(2.0);
    }
    catch (...) {
        return std::numeric_limits<double>::infinity();
    }
}


double sci3545_1(double x) {
    double z = 1.0 / sqrt(x);
    double cdf = gsl_cdf_ugaussian_P(z);  // Φ(z)
    return 2.0 * (1.0 - cdf);
}

double ei(double x) {
    // double val = sin(boost::math::erf_inv(x))/(cos(x)-exp(x));
    double eix = gsl_cdf_ugaussian_Pinv(x/2.0 + 0.5) / sqrt(2.0);
    double val = sin(eix) / (cos(x) - exp(x));
    return val;
}

double Q1_W(double x) {
    double lx = gsl_sf_lambert_W0(x);
    double val = (1+gsl_sf_legendre_Q1(x))/(lx*lx);
    return val;
}

double bj_tan(double x) {
    double jx = gsl_sf_bessel_j0(x);
    double val = (1-jx)/(x*tan(x));
    return val;
}

double Si_tan(double x) {
    double sx = gsl_sf_Si(x);
    double val = (sx - tan(x)) / (x * x * x);
    return val;
}

double by_psi(double x) {
    double yx = (x >= 0) ? gsl_sf_bessel_y0(x) : -gsl_sf_bessel_y0(-x);
    double tx = gsl_sf_psi_1(x);
    double val = yx * yx - tx;
    return val;
}

double fdm_log(double x) {
    double fdx = gsl_sf_fermi_dirac_m1(x);
    double val = (2*fdx-1)/log(1+x);
    return val;
}

double eQ_sqrt(double x) {
    double qx = gsl_sf_erf_Q(x);
    double val = (2 * qx - sqrt(1 + x)) / x;
    return val;
}

double W_var(double x) {
    double lx = gsl_sf_lambert_W0(x);
    return (lx-1) / (lx*lx - 1);
}

double W_log(double x) {
    double lx = gsl_sf_lambert_W0(x);
    return (lx-1) / (lx*log(x) - 1);
}

double pow_df(double x) {
    double val = pow(1+alglib::dawsonintegral(x), 1.0/x);
    return val;
}

double chi_ci(double x) {
    double chix, shix;
    alglib::hyperbolicsinecosineintegrals(x, shix, chix);
    double cix, six;
    alglib::sinecosineintegrals(x, six, cix);
    double val = (chix - cix)/(x*x);
    return val;
}

double fc_bj(double x) {
    double sx, cx;
    alglib::fresnelintegral(x, cx, sx);
    double val = 1.0/cx + gsl_sf_bessel_j1(x) - sin(x)/(x*x);
    return val;
}

double gb_sqrt(double x) {
    double v1 = gsl_sf_gegenpoly_1(0.5, x);
    double v2 = gsl_sf_gegenpoly_n(0, 1, x);
    double val = (sqrt(v1 + 1) - 2) / (x - 3*v2);
    return val;
}

double l1_l2(double x) {
    double l1x = gsl_sf_laguerre_1(0, x);
    double l2x = gsl_sf_laguerre_2(-1, x);
    double val = (log(x)+l1x) / (l2x+0.5);
    return val;
}

double hyp_g2(double x) {
    double hx = gsl_sf_hypot(x, 1.0/x);
    double gx = gsl_sf_gegenpoly_2(1, x);
    double val = (sqrt(2) * hx - 2.0) / (gx + 5.0 - 8.0*x);
    return val;
}

double cos_x2_ofp_fix(double x) {
    static const double c[] = {1.0/2.0, -1.0/24.0, 1.0/720.0, -1.0/40320.0, 1.0/3628800.0, -1.0/479001600.0};
    return even_degree_series(x, c, 6);
}

double exp_2_ofp_fix(double x) {
    static const double c[] = {1.0, 1.0/12.0, 1.0/360.0, 1.0/20160.0, 1.0/1814400.0};
    return odd_degree_series(x, c, 5);
}

double cos_sin_ofp_fix(double x) {
    static const double c[] = {1.0/2.0, 1.0/24.0, 1.0/240.0, 17.0/40320.0, 31.0/725760.0};
    return odd_degree_series(x, c, 5);
}

double sin_sin_ofp_fix(double x) {
    double eps = 1e-6;
    double terms[] = {
        eps * cos(x),
        -0.5 * pow(eps, 2) * sin(x),
        -(1.0/6.0) * pow(eps, 3) * cos(x),
        (1.0/24.0) * pow(eps, 4) * sin(x),
        (1.0/120.0) * pow(eps, 5) * cos(x),
        -(1.0/720.0) * pow(eps, 6) * sin(x),
        -(1.0/5040.0) * pow(eps, 7) * cos(x),
        (1.0/40320.0) * pow(eps, 8) * sin(x),
        (1.0/362880.0) * pow(eps, 9) * cos(x),
        -(1.0/3628800.0) * pow(eps, 10) * sin(x),
        -(1.0/39916800.0) * pow(eps, 11) * cos(x)
    };
    double sum = 0;
    for (int i = 0; i < expansion_term_count(1, 1, 11); i++) {
        sum += terms[i];
    }
    return sum;
}

double tan_tan_ofp_fix(double x) {
    double eps = 1e-6;
    double terms[] = {
        eps * pow(1 / cos(x), 2),
        pow(eps, 2) * tan(x) * pow(1 / cos(x), 2),
        -(1.0 / 3) * pow(eps, 3) * ((cos(2 * x) - 2) * pow(1 / cos(x), 4)),
        -(1.0 / 6) * pow(eps, 4) * ((cos(2 * x) - 5) * tan(x) * pow(1 / cos(x), 4)),
        (1.0 / 60) * pow(eps, 5) * (-26 * cos(2 * x) + cos(4 * x) + 33) * pow(1 / cos(x), 6),
        (1.0 / 360) * pow(eps, 6) * (302 * sin(x) - 57 * sin(3 * x) + sin(5 * x)) * pow(1 / cos(x), 7),
        -(1.0 / 2520) * pow(eps, 7) * (1191 * cos(2 * x) - 120 * cos(4 * x) + cos(6 * x) - 1208) * pow(1 / cos(x), 8),
        pow(eps, 9) * (-88234 * cos(2 * x) + 14608 * cos(4 * x) - 502 * cos(6 * x) + cos(8 * x) + 78095) * pow(1 / cos(x), 10) / 181440.0
    };
    double sum = 0;
    for (int i = 0; i < 8; i++) {
        static const int degrees[] = {1, 2, 3, 4, 5, 6, 7, 9};
        if (degrees[i] > gExpansionTerms) continue;
        sum += terms[i];
    }
    return sum;
}

double cos_cos_ofp_fix(double x) {
    double eps = 1e-6;
    double terms[] = {
        -eps * sin(x),
        -0.5 * pow(eps, 2) * cos(x),
        (1.0/6) * pow(eps, 3) * sin(x),
        (1.0/24) * pow(eps, 4) * cos(x),
        -(1.0/120) * pow(eps, 5) * sin(x),
        -(1.0/720) * pow(eps, 6) * cos(x),
        (1.0/5040) * pow(eps, 7) * sin(x),
        (1.0/40320) * pow(eps, 8) * cos(x),
        -(1.0/362880) * pow(eps, 9) * sin(x),
        -(1.0/3628800) * pow(eps, 10) * cos(x)
    };
    double sum = 0;
    for (int i = 0; i < expansion_term_count(1, 1, 10); i++) {
        sum += terms[i];
    }
    return sum;
}

double exp_exp_ofp_fix(double x) {
    static const double c[] = {2.0, 1.0/3.0, 1.0/60.0, 1.0/2520.0, 1.0/181440.0, 1.0/19958400.0};
    return odd_degree_series(x, c, 6);
}

double exp_1_ofp_fix(double x) {
    static const double c[] = {1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0, 1.0/40320.0, 1.0/362880.0, 1.0/3628800.0};
    return x * horner_series(x, c, expansion_term_count(1, 1, 10));
}

double x_tan_ofp_fix(double x) {
    static const double c[] = {1.0/3.0, 1.0/45.0, 2.0/945.0, 1.0/4725.0, 2.0/93555.0, 1382.0/638512875.0};
    return odd_degree_series(x, c, 6);
}

double log_log_ofp_fix(double x) {
    static const double c[] = {-1.0, -1.0, -1.0/2.0, -5.0/12.0, -7.0/24.0, -191.0/720.0, -33.0/160.0, -11779.0/60480.0, -19309.0/120960.0, -79771.0/518400.0, -945251.0/7257600.0};
    return ordinary_series(x, c, 11);
}

double log_x_ofp_fix(double x) {
    static const double c[] = {-2.0, -2.0/3.0, -2.0/5.0, -2.0/7.0, -2.0/9.0};
    return odd_degree_series(x, c, 5);
}

double sqrt_exp_ofp_fix(double x) {
    const double r2 = sqrt(2.0);
    const double c[] = {r2, 1.0/(2.0*r2), 3.0/(16.0*r2), 7.0/(192.0*r2), 3.0/(1024.0*r2), 1.0/(61440.0*r2), 41.0/(491520.0*r2), 967.0/(41287680.0*r2), -1637.0/(440401920.0*r2), -9737.0/(6794772480.0*r2), 20809.0/(70464307200.0*r2)};
    return ordinary_series(x, c, 11);
}

double sin_tan_ofp_fix(double x) {
    static const double c[] = {-0.5, 9.0/40.0, -27.0/2800.0, 27.0/112000.0, -201.0/86240000.0, 129.0/2038400000.0};
    return even_degree_series(x, c, 6);
}

double exp_x_ofp_fix(double x) {
    static const double c[] = {1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0, 1.0/40320.0, 1.0/362880.0, 1.0/3628800.0, 1.0/39916800.0};
    return ordinary_series(x, c, 11);
}

double x_x2_ofp_fix(double x) {
    static const double c[] = {1.0/2.0, -1.0/4.0, 1.0/8.0, -1.0/16.0, 1.0/32.0, -1.0/64.0, 1.0/128.0, -1.0/256.0, 1.0/512.0, -1.0/1024.0, 1.0/2048.0};
    return shifted_series(x, 1.0, c, 11);
}

double exp_BI_ofp_fix(double x) {
    static const double c[] = {1.0/2.0, 1.0/8.0, 1.0/12.0, 13.0/384.0, 41.0/3840.0, 211.0/46080.0, 109.0/71680.0, 5209.0/10321920.0, 33139.0/185794560.0, 203641.0/3715891200.0, 725819.0/40874803200.0};
    return ordinary_series(x, c, 11);
}

double bJ_sin_ofp_fix(double x) {
    static const double c[] = {1.0/4.0, 5.0/192.0, 31.0/11520.0, 4247.0/15482880.0, 25861.0/928972800.0};
    return odd_degree_series(x, c, 5);
}

double di_tan_ofp_fix(double x) {
    static const double c[] = {-0.25, 41.0/144.0, -13.0/576.0, 4609.0/518400.0, -6151.0/691200.0, -3135683.0/731566080.0, -70921157.0/14631321600.0, -47354665271.0/13168189440000.0, -162045746237.0/52672757760000.0, -578143816489471.0/229442532802560000.0, -653703948745543.0/305923377070080000.0};
    return ordinary_series(x, c, 11);
}



double log_erf_ofp_fix(double x) {
    double sqrt_pi = std::sqrt(M_PI);
    double c[] = {
        -2.0 / sqrt_pi,
        -(2.0 + sqrt_pi) / M_PI,
        (-16.0 - 6.0 * sqrt_pi + 5.0 * M_PI) / (6.0 * std::pow(M_PI, 1.5)),
        (-48.0 - 16.0 * sqrt_pi + 18.0 * M_PI + 3.0 * std::pow(M_PI, 1.5)) / (12.0 * M_PI * M_PI),
        (-2304.0 - 720.0 * sqrt_pi + 1040.0 * M_PI + 210.0 * std::pow(M_PI, 1.5) - 73.0 * M_PI * M_PI) / (360.0 * std::pow(M_PI, 2.5)),
        -(7680.0 + 2304.0 * sqrt_pi - 4080.0 * M_PI - 880.0 * std::pow(M_PI, 1.5) + 490.0 * M_PI * M_PI + 79.0 * std::pow(M_PI, 2.5)) / (720.0 * std::pow(M_PI, 3.0)),
        (-552960.0 - 161280.0 * sqrt_pi + 338688.0 * M_PI + 75600.0 * std::pow(M_PI, 1.5) - 55664.0 * M_PI * M_PI - 8862.0 * std::pow(M_PI, 2.5) + 2275.0 * std::pow(M_PI, 3.0)) / (30240.0 * std::pow(M_PI, 3.5)),
        (-1935360.0 - 552960.0 * sqrt_pi + 1344000.0 * M_PI + 306432.0 * std::pow(M_PI, 1.5) - 278544.0 * M_PI * M_PI - 47376.0 * std::pow(M_PI, 2.5) + 16558.0 * std::pow(M_PI, 3.0) + 317.0 * std::pow(M_PI, 3.5)) / (60480.0 * std::pow(M_PI, 4.0)),
        -512.0 / (9.0 * std::pow(M_PI, 4.5)) - 16.0 / std::pow(M_PI, 4.0) + 928.0 / (21.0 * std::pow(M_PI, 3.5)) + 92.0 / (9.0 * std::pow(M_PI, 3.0)) - 2482.0 / (225.0 * std::pow(M_PI, 2.5)) - 143.0 / (72.0 * M_PI * M_PI) + 4223.0 / (4536.0 * std::pow(M_PI, 1.5)) + 29.0 / (320.0 * M_PI) + 2269.0 / (1814400.0 * sqrt_pi),
        -(371589120.0 + 103219200.0 * sqrt_pi - 319334400.0 * M_PI - 74649600.0 * std::pow(M_PI, 1.5) + 93166080.0 * M_PI * M_PI + 17531136.0 * std::pow(M_PI, 2.5) - 10325520.0 * std::pow(M_PI, 3.0) - 1308560.0 * std::pow(M_PI, 3.5) + 265646.0 * std::pow(M_PI, 4.0) + 52989.0 * std::pow(M_PI, 4.5)) / (3628800.0 * std::pow(M_PI, 5.0)),
        -2048.0 / (11.0 * std::pow(M_PI, 5.5)) - 256.0 / (5.0 * std::pow(M_PI, 5.0)) + 4736.0 / (27.0 * std::pow(M_PI, 4.5)) + 124.0 / (3.0 * std::pow(M_PI, 4.0)) - 3688.0 / (63.0 * std::pow(M_PI, 3.5)) - 57.0 / (5.0 * std::pow(M_PI, 3.0)) + 76243.0 / (9450.0 * std::pow(M_PI, 2.5)) + 7037.0 / (6048.0 * M_PI * M_PI) - 503483.0 / (1360800.0 * std::pow(M_PI, 1.5)) - 27509.0 / (725760.0 * M_PI) + 244247.0 / (21772800.0 * sqrt_pi)
    };
    return ordinary_series(x, c, 11);
}


double acos_fd_ofp_fix(double x) {
    static const double c[] = {-M_PI - std::log(8.0), 1.0/4.0, -1.0/8.0 - M_PI/6.0, 1.0/3.0, 1.0/320.0 - 3.0*M_PI/40.0, 8.0/45.0, -1.0/6720.0 - 5.0*M_PI/112.0, 4.0/35.0, 17.0/1935360.0 - 35.0*M_PI/1152.0, 128.0/1575.0, -31.0/53222400.0 - 63.0*M_PI/2816.0};
    return ordinary_series(x, c, 11);
}

double sci3545_2_ofp_fix(double q) {
    double p = q / 2.0;

    int n = ofp_series_terms(gExpansionTerms, 4);
    double L = -log(p);
    double w = sqrt(2.0 * L);
    double A = log(w) + 0.5 * log(2.0 * M_PI);
    double val = w;

    if (n >= 2) {
        val += -A / w;
    }
    if (n >= 3) {
        val += (-0.5 * A * A + A - 1.0) / (w * w * w);
    }
    if (n >= 4) {
        val += (-0.5 * A * A * A + 2.0 * A * A - 4.0 * A + 3.5)
               / (w * w * w * w * w);
    }

    return val;
}

double sci4034_1_ofp_fix(double x) {
    static const double c[] = {0.5, -0.125, 1.0/48.0, -1.0/384.0, 1.0/3840.0, -1.0/46080.0};
    return x * x * even_series(x, c, expansion_term_count(2, 2, 6));
}

double sci4034_2_ofp_fix(double x) {
    static const double c[] = {2.0, 1.0, 2.0/3.0, 0.5, 0.4, 1.0/3.0, 2.0/7.0, 0.25, 2.0/9.0, 0.2};
    return sqrt(x * horner_series(x, c, expansion_term_count(0, 1, 10)));
}

double sci3547_ofp_fix(double y) {
    int terms = ofp_series_terms(gExpansionTerms, 4);
    return detail::erfc_inv_small_y(y, terms);
}

double sci3545_1_ofp_fix(double x) {
    double z = sqrt(0.5 / x);  
    double sqrt_pi = sqrt(M_PI);
    double z2 = z * z;

    double sum = 1.0;
    double term = 1.0;

    for (int k = 1; k < expansion_term_count(0, 1, 20); ++k) {
        term *= -(2.0 * k - 1.0) / (2.0 * z2);
        sum += term;
    }

    return (exp(-z2) / (z * sqrt_pi)) * sum;
}

double ei_ofp_fix(double x) {
    double sqrt_pi = std::sqrt(M_PI);
    double c[] = {
        -sqrt_pi / 2.0,
        sqrt_pi / 2.0,
        -sqrt_pi * (20.0 + M_PI) / 48.0,
        sqrt_pi * (16.0 + M_PI) / 48.0,
        -sqrt_pi * (2992.0 + 200.0 * M_PI + 27.0 * std::pow(M_PI, 2.0)) / 11520.0,
        sqrt_pi * (2320.0 + 160.0 * M_PI + 27.0 * std::pow(M_PI, 2.0)) / 11520.0,
        -sqrt_pi * (301760.0 + 20944.0 * M_PI + 3780.0 * std::pow(M_PI, 2.0) + 651.0 * std::pow(M_PI, 3.0)) / 1935360.0,
        sqrt_pi * (233472.0 + 16240.0 * M_PI + 3024.0 * std::pow(M_PI, 2.0) + 651.0 * std::pow(M_PI, 3.0)) / 1935360.0,
        -sqrt_pi * (86697216.0 + 6035200.0 * M_PI + 1130976.0 * std::pow(M_PI, 2.0) + 260400.0 * std::pow(M_PI, 3.0) + 47925.0 * std::pow(M_PI, 4.0)) / 928972800.0,
        sqrt_pi * (13413632.0 + 933888.0 * M_PI + 175392.0 * std::pow(M_PI, 2.0) + 41664.0 * std::pow(M_PI, 3.0) + 9585.0 * std::pow(M_PI, 4.0)) / 185794560.0,
        -sqrt_pi * (13696998400.0 + 953669376.0 * M_PI + 179245440.0 * std::pow(M_PI, 2.0) + 42851424.0 * std::pow(M_PI, 3.0) + 10543500.0 * std::pow(M_PI, 4.0) + 1831731.0 * std::pow(M_PI, 5.0)) / 245248819200.0
    };
    return ordinary_series(x, c, 11);
}

double Q1_W_ofp_fix(double x) {
    static const double c[] = {1.0, 2.0, 1.0/3.0, 1.0, -7.0/15.0, 67.0/36.0, -307.0/105.0, 2521.0/360.0, -14039.0/945.0, 31188749.0/907200.0, -493817.0/6237.0};
    return ordinary_series(x, c, 11);
}

double bj_tan_ofp_fix(double x) {
    static const double c[] = {1.0/6.0, -23.0/360.0, -11.0/15120.0, -143.0/604800.0, -361.0/17107200.0, -1416533.0/653837184000.0};
    return even_degree_series(x, c, 6);
}

double Si_tan_ofp_fix(double x) {
    static const double c[] = {-7.0/18.0, -79.0/600.0, -127.0/2352.0, -71423.0/3265920.0, -555959.0/62726400.0, -3589967.0/999398400.0};
    return even_degree_series(x, c, 6);
}

double by_psi_ofp_fix(double x) {
    double pi = M_PI;
    double pi2 = pi * pi;
    double pi4 = pi2 * pi2;
    double pi6 = pi4 * pi2;
    double pi8 = pi6 * pi2;
    double pi10 = pi8 * pi2;
    double pi12 = pi10 * pi2;
    double c[] = {
        -1.0 - pi2 / 6.0,
        -polygamma_even_at_one(2),
        (1.0 / 3.0) - pi4 / 30.0,
        -polygamma_even_at_one(4) / 6.0,
        (-42.0 - 5.0 * pi6) / 945.0,
        -polygamma_even_at_one(6) / 120.0,
        (30.0 - 7.0 * pi8) / 9450.0,
        -polygamma_even_at_one(8) / 5040.0,
        (-22.0 - 15.0 * pi10) / 155925.0,
        -polygamma_even_at_one(10) / 362880.0,
        (2730.0 - 7601.0 * pi12) / 638512875.0
    };
    return ordinary_series(x, c, 11);
}

double fdm_log_ofp_fix(double x) {
    static const double c[] = {1.0/2.0, 1.0/4.0, -1.0/12.0, 0.0, -1.0/180.0, 7.0/720.0, -823.0/120960.0, 1177.0/241920.0, -3319.0/806400.0, 10319.0/2903040.0, -65129.0/21288960.0};
    return ordinary_series(x, c, 11);
}

double eQ_sqrt_ofp_fix(double x) {
    double sqrt_2pi = std::sqrt(2.0 * M_PI);
    double c[] = {
        -std::sqrt(2.0 / M_PI) - 0.5,
        1.0/8.0,
        -1.0/16.0 + 1.0/(3.0 * sqrt_2pi),
        5.0/128.0,
        -1.0/(20.0 * sqrt_2pi) - 7.0/256.0,
        21.0/1024.0,
        -33.0/2048.0 + 1.0/(168.0 * sqrt_2pi),
        429.0/32768.0,
        -1.0/(1728.0 * sqrt_2pi) - 715.0/65536.0,
        2431.0/262144.0,
        -4199.0/524288.0 - 1.0/(21120.0 * sqrt_2pi)
    };
    return ordinary_series(x, c, 11);
}

double W_var_ofp_fix(double x) {
    double delta = x - M_E;
    double inv_e = 1.0 / M_E;
    double c[] = {0.5, -inv_e/8.0, 5.0*inv_e*inv_e/64.0, -43.0*std::pow(inv_e, 3.0)/768.0, 523.0*std::pow(inv_e, 4.0)/12288.0, -2739.0*std::pow(inv_e, 5.0)/81920.0, 158197.0*std::pow(inv_e, 6.0)/5898240.0, -3605171.0*std::pow(inv_e, 7.0)/165150720.0};
    return horner_series(delta, c, expansion_term_count(0, 1, 8));
}

double W_log_ofp_fix(double x) {
    const double e = exp(1.0);
    double delta = x - e;
    double c[] = {0.3333333333333333, -0.030656620097620196, 0.007675264905780119, -0.002204115005868977, 0.0006707681714134265, -0.00021077950097618563, 6.759635969842329e-05, -2.198607515321714e-05, 7.225880712334929e-06, -2.39397506957339e-06, 1e-6};
    return horner_series(delta, c, expansion_term_count(0, 1, 11));
}

double pow_df_ofp_fix(double x) {
    double e = exp(1.0);
    double c[] = {1.0, -1.0/2.0, -5.0/24.0, 9.0/16.0, -2257.0/5760.0, 37.0/11520.0, 24445.0/82944.0, -1989077.0/5806080.0, 236533639.0/1393459200.0, 22039481.0/309657600.0, -1489734851.0/6688604160.0};
    return e * ordinary_series(x, c, 11);
}

double chi_ci_ofp_fix(double x) {
    static const double c[] = {0.5, 1.0/2160.0, 1.0/18144000.0, 1.0/610248038400.0};
    double x4 = x * x * x * x;
    return horner_series(x4, c, expansion_term_count(0, 4, 4));
}

double fc_bj_ofp_fix(double x) {
    double c[] = {1.0/2.0, (3.0*M_PI*M_PI - 5.0)/120.0, 1.0/720.0, (203.0*std::pow(M_PI, 4.0) - 15.0)/604800.0, 1.0/3628800.0};
    return odd_degree_series(x, c, 5);
}

double sympy_28514_log_ratio(double b) {
    return std::log((b / (b + 1.0)) / b);
}

double sympy_12671_log_cosh(double x) {
    return std::log(std::cosh(x));
}

double statsmodels_1604_logistic(double eta) {
    double ex = std::exp(eta);
    return ex / (1.0 + ex);
}

double statsmodels_826_cloglog(double p) {
    return std::log(-std::log(1.0 - p));
}

double scipy_17194_expm1_kernel(double g) {
    return 1.0 - std::exp(-g);
}

double pytorch_39242_log1mexp(double x) {
    return std::log(1.0 - std::exp(x));
}

double sympy_28514_log_ratio_ofp_fix(double b) {
    return -log1p_expansion(b, 16);
}

double sympy_12671_log_cosh_ofp_fix(double x) {
    double ax = std::fabs(x);
    return ax - std::log(2.0) + log1p_expansion(std::exp(-2.0 * ax), 16);
}

double statsmodels_1604_logistic_ofp_fix(double eta) {
    double q = std::exp(-eta);
    double total = 1.0;
    double power = 1.0;
    int terms = expansion_term_count(1, 1, 12);
    for (int n = 1; n <= terms; ++n) {
        power *= q;
        total += (n % 2 ? -1.0 : 1.0) * power;
    }
    return total;
}

double statsmodels_826_cloglog_ofp_fix(double p) {
    double inner = 0.0;
    double power = 1.0;
    int terms = expansion_term_count(1, 1, 16);
    for (int n = 1; n <= terms; ++n) {
        power *= p;
        inner += power / n;
    }
    return std::log(inner);
}

double scipy_17194_expm1_kernel_ofp_fix(double g) {
    double total = 0.0;
    double power = 1.0;
    double factorial = 1.0;
    int terms = expansion_term_count(1, 1, 18);
    for (int n = 1; n <= terms; ++n) {
        power *= g;
        factorial *= n;
        total += (n % 2 ? 1.0 : -1.0) * power / factorial;
    }
    return total;
}

double pytorch_39242_log1mexp_ofp_fix(double x) {
    return std::log(-expm1_expansion(x, 24));
}

