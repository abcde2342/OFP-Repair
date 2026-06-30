#include <cmath>
#include <mpreal.h>
#include <vector>
#include <random>
#include <gsl/gsl_sf.h>
#include <cstdio>
#include <functional>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>


using mpfr::mpreal;
// using namespace std;

// typedef long double real;
typedef mpreal real;

const mpfr_prec_t highprec=1024;

real naive_trans_gsl_lambert_W0(real x);
real naive_trans_gsl_sf_legendre_Q1(real x);

namespace {
///////////////////////////////////////////////////////////////////////////////
// GSL Utility
///////////////////////////////////////////////////////////////////////////////
// #define INT_MIN -2147483648

#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)
// #define M_E        2.71828182845904523536028747135      /* e */
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */

#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */

#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

#define GSL_SUCCESS 0
#define GSL_FAILURE 1

#define GSL_NAN NAN

#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

real
GSL_MAX_DBL (real a, real b)
{
  return GSL_MAX (a, b);
}
real
GSL_MIN_DBL (real a, real b)
{
  return GSL_MIN (a, b);
}

struct gsl_sf_result_struct {
  real val;
  real err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   return result.val;


#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

///////////////////////////////////////////////////////////////////////////////
// GSL Lambert W Function
///////////////////////////////////////////////////////////////////////////////

static real halley_iteration(real x, real w_initial, unsigned int max_iters) {
    real w = w_initial;
    unsigned int i;

    for (i = 0; i < max_iters; i++) {
        real tol;
        const real e = exp(w);
        const real p = w + 1.0;
        real t = w*e - x;

        if (w > 0) {
            t = (t/p)/e;
        } else {
            t /= e*p - 0.5*(p + 1.0)*t/p;
        };

        w -= t;

        tol = 10 * GSL_DBL_EPSILON * GSL_MAX_DBL(fabs(w), 1.0/(fabs(p)*e));

        if (fabs(t) < tol) {
            return w;
        }
    }
    // Should never get here
    return w;
}

static real
series_eval(real r)
{
  static const real c[12] = {
    -1.0,
     2.331643981597124203363536062168,
    -1.812187885639363490240191647568,
     1.936631114492359755363277457668,
    -2.353551201881614516821543561516,
     3.066858901050631912893148922704,
    -4.175335600258177138854984177460,
     5.858023729874774148815053846119,
    -8.401032217523977370984161688514,
     12.250753501314460424,
    -18.100697012472442755,
     27.029044799010561650
  };
  const real t_8 = c[8] + r*(c[9] + r*(c[10] + r*c[11]));
  const real t_5 = c[5] + r*(c[6] + r*(c[7]  + r*t_8));
  const real t_1 = c[1] + r*(c[2] + r*(c[3]  + r*(c[4] + r*t_5)));
  return c[0] + r*t_1;
}
} // namespace

real naive_trans_gsl_lambert_W0(real x) {
    const real one_over_E = 1.0/M_E;
    const real q = x + one_over_E;

    if (x == 0.0) {
        return 0.0;
    }
    else if (q < 0.0) {
        return -1.0;
    }
    else if (q == 0.0) {
        return -1.0;
    }
    else if (q < 1.0e-3) {
        const real r = sqrt(q);
        real val = series_eval(r);
        return val;
    }
    else {
        static const unsigned int MAX_ITERS = 10;
        real w;

        if (x < 1.0) {
            const real p = sqrt(2.0 * M_E * q);
            w = -1.0 + p*(1.0 + p*(-1.0/3.0 + p*11.0/72.0));
        }
        else {
            w = log(x);
            if (x > 3.0) w -= log(w);
        }

        return halley_iteration(x, w, MAX_ITERS);
    }
}

real naive_trans_gsl_sf_legendre_Q1(real x)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0 || x == 1.0) {
    return 0;
  }
  else if(x*x < GSL_ROOT6_DBL_EPSILON) { /* |x| <~ 0.05 */
    const real c3 = 1.0/3.0;
    const real c5 = 1.0/5.0;
    const real c7 = 1.0/7.0;
    const real c9 = 1.0/9.0;
    const real c11 = 1.0/11.0;
    const real y = x * x;
    const real series = 1.0 + y*(c3 + y*(c5 + y*(c7 + y*(c9 + y*c11))));
    return x * x * series - 1.0;
  }
  else if(x < 1.0){
    return 0.5 * x * (log((1.0+x)/(1.0-x))) - 1.0;
  }
  else if(x < 6.0) {
    return 0.5 * x * log((x+1.0)/(x-1.0)) - 1.0;
  }
  else if(x*GSL_SQRT_DBL_MIN < 0.99/M_SQRT3) {
    const real y = 1/(x*x);
    const real c1 = 3.0/5.0;
    const real c2 = 3.0/7.0;
    const real c3 = 3.0/9.0;
    const real c4 = 3.0/11.0;
    const real c5 = 3.0/13.0;
    const real c6 = 3.0/15.0;
    const real c7 = 3.0/17.0;
    const real c8 = 3.0/19.0;
    const real sum = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*c8)))))));
    return sum / (3.0*x*x);
  }
  else {
    return 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
// First Example: fun 26 with fake lifted oracle
///////////////////////////////////////////////////////////////////////////////
double fun26(double x) {
    double lx = gsl_sf_lambert_W0(x);
    double val = (1+gsl_sf_legendre_Q1(x))/(lx*lx);
    return val;
}

real fake_oracle_26(real x) {
    if (std::is_same<real, mpfr::mpreal>::value) {
        mpfr::mpreal::set_default_prec(highprec);
    }
    real lx = naive_trans_gsl_lambert_W0(x);
    real val = (1+naive_trans_gsl_sf_legendre_Q1(x))/(lx*lx);
    return val;
}

mpfr::mpreal oracle26(double dx) {
    // https://www.wolframalpha.com/input?i=series+%281%2BLegendreQ%5B1%2Cx%5D%29%2F%28W%28x%29%5E2%29+at+x+%3D+0
    // https://math.stackexchange.com/questions/1700919/how-to-derive-the-lambert-w-function-series-expansion
    mpfr::mpreal::set_default_prec(highprec);
    mpfr::mpreal x = dx, y;
    int num = 12;
    mpfr::mpreal c[12] = {
        1.0,
        2.0,
        ((mpfr::mpreal)1.0) / 3.0,
        1.0,
        -((mpfr::mpreal)7.0) / 15.0,
        ((mpfr::mpreal)67.0) / 36.0,
        -((mpfr::mpreal)307.0) / 105.0,
        ((mpfr::mpreal)2521.0) / 360.0,
        -((mpfr::mpreal)14039.0) / 945.0,
        mpfr::mpreal("31188749") / mpfr::mpreal("907200"),
        -mpfr::mpreal("493817") / mpfr::mpreal("6237"),
        mpfr::mpreal("1017580243") / mpfr::mpreal("5443200"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

///////////////////////////////////////////////////////////////////////////////
// Second Example: variate_lambert with fake lifted oracle
///////////////////////////////////////////////////////////////////////////////
double variate_lambert(double x) {
    double lx = gsl_sf_lambert_W0(x);
    return (lx-1) / (lx*lx - 1);
}

real fake_oracle_variate_lambert(real x) {
    real lx = naive_trans_gsl_lambert_W0(x);
    return (lx-1) / (lx*lx - 1);
}

mpreal oracle_variate_lambert(double dx) {
    boost::multiprecision::cpp_dec_float_50 x = dx, lx, y;
    lx = boost::math::lambert_w0(x);
    y = 1 / (lx + 1);
    return (double) y;
}

///////////////////////////////////////////////////////////////////////////////
// Third Example: variate_lambert_2 with fake lifted oracle
///////////////////////////////////////////////////////////////////////////////
double variate_lambert_2(double x) {
    double lx = gsl_sf_lambert_W0(x);
    return (lx-1) / (lx*log(x) - 1);
}

real fake_oracle_variate_lambert_2(real x) {
    real lx = naive_trans_gsl_lambert_W0(x);
    return (lx-1) / (lx*log(x) - 1);
}

mpreal oracle_variate_lambert_2(double dx) {
    boost::multiprecision::cpp_dec_float_50 x = dx, lx, y;
    lx = boost::math::lambert_w0(x);
    y = (lx-1) / (lx*log(x) - 1);
    return (double) y;
}

// std::function<double(double)> originTarget = fun26;
// std::function<real(real)> liftedTarget = fake_oracle_26;
// std::function<mpreal(double)> oracle = oracle26;
// double lbound = -1e-2, rbound = 1e-2;

// std::function<double(double)> originTarget = variate_lambert;
// std::function<real(real)> liftedTarget = fake_oracle_variate_lambert;
// std::function<mpreal(double)> oracle = oracle_variate_lambert;
// double E = 2.71828182845904523536;
// double lbound = E-1e-2, rbound = E+1e-2;

std::function<double(double)> originTarget = variate_lambert_2;
std::function<real(real)> liftedTarget = fake_oracle_variate_lambert_2;
std::function<mpreal(double)> oracle = oracle_variate_lambert_2;
double E = 2.71828182845904523536;
double lbound = E-1e-2, rbound = E+1e-2;

std::vector<double> randomOnValue(double l, double r, int nums) {
    std::vector<double> xs;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937_64 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(l, r);

    for (int i = 0; i < nums; i++) {
        xs.push_back(dis(gen));
    }
    return xs;
}

std::vector<double> randomOnMagnitude(double lowMag, double HighMag, int nums) {
    std::vector<double> xs;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937_64 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(lowMag, HighMag);

    for (int i = 0; i < nums; i++) {
        double x = exp10(dis(gen));
        xs.push_back(x);
    }
    return xs;
}

int evaluationNum = 5000;

void measureErr() {
    // On accurate area
    double stable_origin_max_relative_err = 0;
    double stable_fake_oracle_max_relative_err = 0;

    double stable_origin_max_absolute_err = 0;
    double stable_fake_oracle_max_absolute_err = 0;

    double diff_origin_fake_oracle_max_relative_err = 0;
    double diff_origin_fake_oracle_max_absolute_err = 0;

    auto inputs = randomOnValue(lbound, rbound, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = inputs[i];
        double origin_res = originTarget(xi);
        double fake_oracle_res = (double)liftedTarget(xi);
        mpfr::mpreal oracle_res = oracle(xi);

        // printf("=============\nx: %.3e\n", xi);
        // printf("origin: %.15e\n", origin_res);
        // printf("poly:   %.15e\n", fake_oracle_res);
        // printf("oracle: %.15e\n", (double)(oracle_res));;
        // getchar();

        double origin_absolute_err = (double)fabs(origin_res - oracle_res);
        double fake_oracle_absolute_err  = (double)fabs(fake_oracle_res - oracle_res);
        double origin_relative_err = (double)fabs((origin_res - oracle_res) / oracle_res);
        double fake_oracle_relative_err  = (double)fabs((fake_oracle_res - oracle_res) / oracle_res);

        stable_origin_max_relative_err = std::max(stable_origin_max_relative_err, origin_relative_err);
        stable_fake_oracle_max_relative_err  = std::max(stable_fake_oracle_max_relative_err, fake_oracle_relative_err);
        stable_origin_max_absolute_err = std::max(stable_origin_max_absolute_err, origin_absolute_err);
        stable_fake_oracle_max_absolute_err  = std::max(stable_fake_oracle_max_absolute_err, fake_oracle_absolute_err);

        double diff_relative_err = (double)fabs((origin_res - fake_oracle_res) / fake_oracle_res);
        double diff_absolute_err = (double)fabs(origin_res - fake_oracle_res);
        
        diff_origin_fake_oracle_max_relative_err = std::max(diff_origin_fake_oracle_max_relative_err, diff_relative_err);
        diff_origin_fake_oracle_max_absolute_err = std::max(diff_origin_fake_oracle_max_absolute_err, diff_absolute_err);
    }
    printf("Diff between origin and fake oracle (abs): %.2e\n", diff_origin_fake_oracle_max_relative_err);
    printf("Diff between origin and fake oracle (rel): %.2e\n\n", diff_origin_fake_oracle_max_relative_err);

    printf("Error on Stable Area for Fake Oracle (abs): %.2e -> %.2e\n", stable_origin_max_absolute_err, stable_fake_oracle_max_absolute_err);
    printf("Error on Stable Area for Fake Oracle (rel): %.2e -> %.2e\n", stable_origin_max_relative_err, stable_fake_oracle_max_relative_err);

    // On decayed area
    double decayed_origin_max_relative_err = 0;
    double decayed_fake_oracle_max_relative_err = 0;

    double decayed_origin_max_absolute_err = 0;
    double decayed_fake_oracle_max_absolute_err = 0;

    double decayed_diff_origin_fake_oracle_max_relative_err = 0;
    double decayed_diff_origin_fake_oracle_max_absolute_err = 0;

    auto mags = randomOnMagnitude(-45, -2, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = (lbound + rbound)/2 + ((i%2==0) ? mags[i] : -mags[i]);
        double origin_res = originTarget(xi);
        double fake_oracle_res = (double)liftedTarget(xi);
        mpfr::mpreal oracle_res = oracle(xi);

        double origin_absolute_err = (double)fabs(origin_res - oracle_res);
        double fake_oracle_absolute_err  = (double)fabs(fake_oracle_res - oracle_res);
        double origin_relative_err = (double)fabs((origin_res - oracle_res) / oracle_res);
        double fake_oracle_relative_err  = (double)fabs((fake_oracle_res - oracle_res) / oracle_res);

        decayed_origin_max_relative_err = std::max(decayed_origin_max_relative_err, origin_relative_err);
        decayed_fake_oracle_max_relative_err  = std::max(decayed_fake_oracle_max_relative_err, fake_oracle_relative_err);
        decayed_origin_max_absolute_err = std::max(decayed_origin_max_absolute_err, origin_absolute_err);
        decayed_fake_oracle_max_absolute_err  = std::max(decayed_fake_oracle_max_absolute_err, fake_oracle_absolute_err);

        // if (fake_oracle_relative_err > 1e-1) {
        //   printf("=============\nx: %.18e, err:%.3e\n", xi, fake_oracle_relative_err);
        //   printf("origin: %.15e\n", origin_res);
        //   printf("fake_oracle:  %.15e\n", fake_oracle_res);
        //   printf("oracle: %.15e\n", (double)(oracle_res));
        // }
    }
    printf("Error on Decayed Area for Fake Oracle (abs): %.2e -> %.2e\n", decayed_origin_max_absolute_err, decayed_fake_oracle_max_absolute_err);
    printf("Error on Decayed Area for Fake Oracle (rel): %.2e -> %.2e\n", decayed_origin_max_relative_err, decayed_fake_oracle_max_relative_err);
}

void testInput() {
    while (1) {
        double d;
        printf("Enter a specific input for manually testing: ");
        scanf("%lf", &d);
        printf("original:           %.15e\n", (double)originTarget(d));
        printf("lifted_fake_oracle: %.15e\n", (double)liftedTarget(d));
        printf("true oracle:        %.15e\n", (double)oracle(d));
    }
}
int main() {
    measureErr();

    testInput();
    return 0;
}