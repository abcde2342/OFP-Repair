#include "gsl-fit.h"
#include <cstdio>
#include "../inputs/oracles.h"
#include <chrono>
#include <random>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>
#include "serializer.h"
#include <array>
#include <limits>
#include <iostream>
#include <sstream>

#define FullInfoLevel 4
#define NoInfoLevel 0

#define USE_WEIGHT 1
#define NO_WEIGHT  0

using namespace std;

extern "C" double autornp_Q1_W(double);
extern "C" double autornp_Si_tan(double);
extern "C" double autornp_W_log(double);
extern "C" double autornp_W_var(double);
extern "C" double autornp_acos_fd(double);
extern "C" double autornp_bj_tan(double);
extern "C" double autornp_by_psi(double);
extern "C" double autornp_chi_ci(double);
extern "C" double autornp_cos_cos(double);
extern "C" double autornp_cos_x2(double);
extern "C" double autornp_di_tan(double);
extern "C" double autornp_eQ_sqrt(double);
extern "C" double autornp_ei(double);
extern "C" double autornp_exp_BI(double);
extern "C" double autornp_fdm_log(double);
extern "C" double autornp_log_erf(double);
extern "C" double autornp_pow_df(double);
extern "C" double autornp_sci3545_2(double);
extern "C" double autornp_sci4034_1(double);
extern "C" double autornp_sci3547(double);
extern "C" double autornp_scipy_17194_expm1_kernel(double);
extern "C" double autornp_sin_sin(double);
extern "C" double autornp_sympy_12671_log_cosh(double);
extern "C" double autornp_tan_tan(double);
extern "C" double autornp_x_x2(double);

extern "C" double herbie_cos_x2(double);
extern "C" double herbie_exp_2(double);
extern "C" double herbie_cos_sin(double);
extern "C" double herbie_sin_sin(double);
extern "C" double herbie_tan_tan(double);
extern "C" double herbie_cos_cos(double);
extern "C" double herbie_exp_exp(double);
extern "C" double herbie_exp_1(double);
extern "C" double herbie_x_tan(double);
extern "C" double herbie_log_log(double);
extern "C" double herbie_log_x(double);
extern "C" double herbie_sqrt_exp(double);
extern "C" double herbie_sin_tan(double);
extern "C" double herbie_exp_x(double);
extern "C" double herbie_x_x2(double);
extern "C" double herbie_sci4034_1(double);
extern "C" double herbie_sci4034_2(double);
extern "C" double herbie_sympy_28514_log_ratio(double);
extern "C" double herbie_sympy_12671_log_cosh(double);
extern "C" double herbie_statsmodels_1604_logistic(double);
extern "C" double herbie_statsmodels_826_cloglog(double);
extern "C" double herbie_scipy_17194_expm1_kernel(double);
extern "C" double herbie_pytorch_39242_log1mexp(double);

bool fileExistsAndNotEmpty(const char* path);

bool isFinitePolynomial(const Polynomial& poly) {
    if (!std::isfinite(poly.lbound) || !std::isfinite(poly.rbound) ||
        !std::isfinite(poly.c1) || !std::isfinite(poly.c2)) {
        return false;
    }
    for (double c : poly.coef) {
        if (!std::isfinite(c)) {
            return false;
        }
    }
    double mid = (poly.lbound + poly.rbound) / 2.0;
    return std::isfinite(poly.evaluate(poly.lbound)) &&
           std::isfinite(poly.evaluate(mid)) &&
           std::isfinite(poly.evaluate(poly.rbound));
}

double safeAbsoluteError(double value, const mpreal& oracle) {
    if (!std::isfinite(value)) {
        return std::numeric_limits<double>::infinity();
    }
    double absOracle = (double)fabs(oracle);
    if (!std::isfinite(absOracle)) {
        return std::numeric_limits<double>::infinity();
    }
    double err = (double)fabs(value - oracle);
    return std::isfinite(err) ? err : std::numeric_limits<double>::infinity();
}

double relativeErrorFloor(int funcIndex) {
    (void)funcIndex;
    return 1e-23;
}

double safeRelativeError(double value, const mpreal& oracle, int funcIndex) {
    if (!std::isfinite(value)) {
        return std::numeric_limits<double>::infinity();
    }
    double absOracle = (double)fabs(oracle);
    if (!std::isfinite(absOracle)) {
        return std::numeric_limits<double>::infinity();
    }
    double denom = std::max(absOracle, relativeErrorFloor(funcIndex));
    double err = (double)fabs(value - oracle) / denom;
    return std::isfinite(err) ? err : std::numeric_limits<double>::infinity();
}

double decayedAbsoluteError(double value, const mpreal& oracle) {
    if (!std::isfinite(value)) {
        return 1.0;
    }
    double absOracle = (double)fabs(oracle);
    if (!std::isfinite(absOracle)) {
        return 1.0;
    }
    double err = (double)fabs(value - oracle);
    return std::isfinite(err) ? err : 1.0;
}

double decayedRelativeError(double value, const mpreal& oracle, int funcIndex) {
    if (!std::isfinite(value)) {
        return 1.0;
    }
    double absOracle = (double)fabs(oracle);
    if (!std::isfinite(absOracle)) {
        return 1.0;
    }
    double denom = std::max(absOracle, relativeErrorFloor(funcIndex));
    double err = (double)fabs(value - oracle) / denom;
    return std::isfinite(err) ? err : 1.0;
}

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

bool hasPositiveLeftEndpointSingularity(int funcIndex) {
    return funcIndex == 32 || funcIndex == 34 || funcIndex == 35 || funcIndex == 36
        || funcIndex == 40;
}

bool hasPositiveRightEndpointSingularity(int funcIndex) {
    return funcIndex == 42;
}

bool shouldSkipAcesoSynthesis(int funcIndex) {
    return funcIndex == 36;
}

const std::vector<std::string> functionIDs = {
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
    "c1 sci3545_2",
    "c2 sci4034_1",
    "c3 sci4034_2",
    "c4 sci3547",
    "c5 sci3545_1",
    "case sympy/sympy#28514",
    "case sympy/sympy#12671",
    "case statsmodels/statsmodels#1604",
    "case statsmodels/statsmodels#826",
    "case scipy/scipy#17194",
    "case pytorch/pytorch#39242",
};

const std::vector<std::pair<double, double>> errIntervals = {
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {2.718281828459045-1e-2, 2.718281828459045+1e-2}, // W_var
    {2.718281828459045-1e-2, 2.718281828459045+1e-2}, // W_log
    {-1e-2, 1e-2}, // 15
    {-1e-2, 1e-2}, // 16
    {-1e-2, 1e-2}, // 17
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {1.5607963267948966192313216916398, 1.5807963267948966192313216916398},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {-1e-2, 1e-2},
    {1-1e-2, 1+1e-2},
    {0, 1e-10},  // c1 sci3545_2
    {-1e-2, 1e-2},  // c2 sci4034_1
    {0, 1e-2},  // c3 sci4034_2
    {0, 1e-10},  // c4 sci3547
    {0, 1e-2},  // c5 sci3545_1
    {0, 1e-2},  // sympy/sympy#28514
    {1000.0 - 1e-2, 1000.0 + 1e-2},  // sympy/sympy#12671
    {1741.90624 - 1e-2, 1741.90624 + 1e-2},  // statsmodels/statsmodels#1604
    {0, 1e-2},  // statsmodels/statsmodels#826
    {-1e-2, 1e-2},  // scipy/scipy#17194
    {-1e-2, 0},  // pytorch/pytorch#39242
};

struct AutoRNPEntry {
    int funcIndex;
    std::function<double(double)> patchFunc;
};

struct HerbieEntry {
    int funcIndex;
    std::function<double(double)> patchFunc;
};

const std::vector<AutoRNPEntry> autoRNPFuncList = {
    {0, autornp_exp_BI},
    {2, autornp_di_tan},
    {3, autornp_log_erf},
    {4, autornp_acos_fd},
    {5, autornp_ei},
    {6, autornp_Q1_W},
    {7, autornp_bj_tan},
    {8, autornp_Si_tan},
    {9, autornp_by_psi},
    {10, autornp_fdm_log},
    {11, autornp_eQ_sqrt},
    {12, autornp_W_var},
    {13, autornp_W_log},
    {14, autornp_pow_df},
    {15, autornp_chi_ci},
    {17, autornp_cos_x2},
    {20, autornp_sin_sin},
    {21, autornp_tan_tan},
    {22, autornp_cos_cos},
    {31, autornp_x_x2},
    {32, autornp_sci3545_2},
    {33, autornp_sci4034_1},
    {35, autornp_sci3547},
    {38, autornp_sympy_12671_log_cosh},
    {41, autornp_scipy_17194_expm1_kernel},
};

const std::vector<HerbieEntry> herbieFuncList = {
    {17, herbie_cos_x2},
    {18, herbie_exp_2},
    {19, herbie_cos_sin},
    {20, herbie_sin_sin},
    {21, herbie_tan_tan},
    {22, herbie_cos_cos},
    {23, herbie_exp_exp},
    {24, herbie_exp_1},
    {25, herbie_x_tan},
    {26, herbie_log_log},
    {27, herbie_log_x},
    {28, herbie_sqrt_exp},
    {29, herbie_sin_tan},
    {30, herbie_exp_x},
    {31, herbie_x_x2},
    {33, herbie_sci4034_1},
    {34, herbie_sci4034_2},
    {37, herbie_sympy_28514_log_ratio},
    {38, herbie_sympy_12671_log_cosh},
    {39, herbie_statsmodels_1604_logistic},
    {40, herbie_statsmodels_826_cloglog},
    {41, herbie_scipy_17194_expm1_kernel},
    {42, herbie_pytorch_39242_log1mexp},
};

const AutoRNPEntry* findAutoRNPEntry(int funcIndex) {
    for (const AutoRNPEntry& entry : autoRNPFuncList) {
        if (entry.funcIndex == funcIndex) {
            return &entry;
        }
    }
    return nullptr;
}

const HerbieEntry* findHerbieEntry(int funcIndex) {
    for (const HerbieEntry& entry : herbieFuncList) {
        if (entry.funcIndex == funcIndex) {
            return &entry;
        }
    }
    return nullptr;
}

const int defaultRootsNum = 1000;
const int defaultDegree   = 6;
const int defaultPoly     = MODE_CHEBY;

const int evaluationNum  = 5000;

struct RepairMetrics {
    double stable_origin_abs = 0;
    double stable_aceso_abs = 0;
    double stable_ofp_repair_abs = 0;
    double stable_origin_rel = 0;
    double stable_aceso_rel = 0;
    double stable_ofp_repair_rel = 0;
    double decayed_origin_abs = 0;
    double decayed_aceso_abs = 0;
    double decayed_ofp_repair_abs = 0;
    double decayed_origin_rel = 0;
    double decayed_aceso_rel = 0;
    double decayed_ofp_repair_rel = 0;

    void add(const RepairMetrics& other) {
        stable_origin_abs += other.stable_origin_abs;
        stable_aceso_abs += other.stable_aceso_abs;
        stable_ofp_repair_abs += other.stable_ofp_repair_abs;
        stable_origin_rel += other.stable_origin_rel;
        stable_aceso_rel += other.stable_aceso_rel;
        stable_ofp_repair_rel += other.stable_ofp_repair_rel;
        decayed_origin_abs += other.decayed_origin_abs;
        decayed_aceso_abs += other.decayed_aceso_abs;
        decayed_ofp_repair_abs += other.decayed_ofp_repair_abs;
        decayed_origin_rel += other.decayed_origin_rel;
        decayed_aceso_rel += other.decayed_aceso_rel;
        decayed_ofp_repair_rel += other.decayed_ofp_repair_rel;
    }

    void divide(double value) {
        stable_origin_abs /= value;
        stable_aceso_abs /= value;
        stable_ofp_repair_abs /= value;
        stable_origin_rel /= value;
        stable_aceso_rel /= value;
        stable_ofp_repair_rel /= value;
        decayed_origin_abs /= value;
        decayed_aceso_abs /= value;
        decayed_ofp_repair_abs /= value;
        decayed_origin_rel /= value;
        decayed_aceso_rel /= value;
        decayed_ofp_repair_rel /= value;
    }
};

struct AutoRNPMetrics {
    double stable_origin_abs = 0;
    double stable_autornp_abs = 0;
    double stable_origin_rel = 0;
    double stable_autornp_rel = 0;
    double decayed_origin_abs = 0;
    double decayed_autornp_abs = 0;
    double decayed_origin_rel = 0;
    double decayed_autornp_rel = 0;

    void add(const AutoRNPMetrics& other) {
        stable_origin_abs += other.stable_origin_abs;
        stable_autornp_abs += other.stable_autornp_abs;
        stable_origin_rel += other.stable_origin_rel;
        stable_autornp_rel += other.stable_autornp_rel;
        decayed_origin_abs += other.decayed_origin_abs;
        decayed_autornp_abs += other.decayed_autornp_abs;
        decayed_origin_rel += other.decayed_origin_rel;
        decayed_autornp_rel += other.decayed_autornp_rel;
    }

    void divide(double value) {
        stable_origin_abs /= value;
        stable_autornp_abs /= value;
        stable_origin_rel /= value;
        stable_autornp_rel /= value;
        decayed_origin_abs /= value;
        decayed_autornp_abs /= value;
        decayed_origin_rel /= value;
        decayed_autornp_rel /= value;
    }
};

int measurePatchRuntime(int funcIndex, int weighted) {
    Polynomial poly;
    double lbound = errIntervals[funcIndex].first;
    double rbound = errIntervals[funcIndex].second;
    int rootsNum  = defaultRootsNum;
    int highestDegree = defaultDegree;
    int polymode = defaultPoly;
    Observation observation(funcIndex, lbound, rbound, rootsNum);
    observation.init();

    GSLFitModule gslfit;
    gslfit.init(highestDegree, weighted, polymode, lbound, rbound);
    gslfit.fitObservations(observation);
    poly = gslfit.getPoly();
    poly.checkacc0();

    // prepare inputs
    int EvalTimeNum = 1000000;
    auto inputsvec = randomOnValue(lbound, rbound, EvalTimeNum);
    double inputs[1000000];
    for (int i = 0; i < EvalTimeNum; i++) {
        inputs[i] = inputsvec[i];
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < EvalTimeNum; i++) {
        poly.evaluate(inputs[i]);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    int64_t time_on_patch = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    auto func = simpleFuncList[funcIndex];
    start_time = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < EvalTimeNum; i++) {
        func(inputs[i]);
    }
    end_time = std::chrono::high_resolution_clock::now();
    int64_t time_on_origin = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    printf("Function and Time: %d, %ld, %ld\n", funcIndex, time_on_patch, time_on_origin);
    return 0;
}

RepairMetrics evaluateRepairMetrics(int funcIndex, const Polynomial& poly, double lbound, double rbound, bool evaluateAceso = true) {
    auto originFunc = simpleFuncList[funcIndex];
    auto fixFunc = fixFuncList[funcIndex];
    auto oracleFunc = oracleFuncList[funcIndex];
    RepairMetrics metrics;

    auto inputs = randomOnValue(lbound, rbound, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = inputs[i];
        double origin_res = originFunc(xi);
        double patch_res = evaluateAceso ? poly.evaluate(xi) : std::numeric_limits<double>::quiet_NaN();
        double ofp_res = fixFunc(xi);
        mpreal oracle_res = oracleFunc(xi);

        double origin_absolute_err = decayedAbsoluteError(origin_res, oracle_res);
        double patch_absolute_err = decayedAbsoluteError(patch_res, oracle_res);
        double ofp_absolute_err = decayedAbsoluteError(ofp_res, oracle_res);
        double origin_relative_err = decayedRelativeError(origin_res, oracle_res, funcIndex);
        double patch_relative_err = decayedRelativeError(patch_res, oracle_res, funcIndex);
        double ofp_relative_err = decayedRelativeError(ofp_res, oracle_res, funcIndex);

        metrics.stable_origin_rel = max(metrics.stable_origin_rel, origin_relative_err);
        metrics.stable_aceso_rel = max(metrics.stable_aceso_rel, patch_relative_err);
        metrics.stable_ofp_repair_rel = max(metrics.stable_ofp_repair_rel, ofp_relative_err);
        metrics.stable_origin_abs = max(metrics.stable_origin_abs, origin_absolute_err);
        metrics.stable_aceso_abs = max(metrics.stable_aceso_abs, patch_absolute_err);
        metrics.stable_ofp_repair_abs = max(metrics.stable_ofp_repair_abs, ofp_absolute_err);
    }

    double decayedHighMag = -2;
    if (hasPositiveLeftEndpointSingularity(funcIndex) || hasPositiveRightEndpointSingularity(funcIndex)) {
        decayedHighMag = log10(rbound - lbound);
    }
    auto mags = randomOnMagnitude(-20, decayedHighMag, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = (lbound + rbound)/2 + ((i%2==0) ? mags[i] : -mags[i]);
        if (hasPositiveLeftEndpointSingularity(funcIndex)) {
            xi = lbound + mags[i];
        } else if (hasPositiveRightEndpointSingularity(funcIndex)) {
            xi = rbound - mags[i];
        }
        double origin_res = originFunc(xi);
        double patch_res = evaluateAceso ? poly.evaluate(xi) : std::numeric_limits<double>::quiet_NaN();
        double ofp_res = fixFunc(xi);
        mpreal oracle_res = oracleFunc(xi);

        double origin_absolute_err = decayedAbsoluteError(origin_res, oracle_res);
        double patch_absolute_err = decayedAbsoluteError(patch_res, oracle_res);
        double ofp_absolute_err = decayedAbsoluteError(ofp_res, oracle_res);
        double origin_relative_err = decayedRelativeError(origin_res, oracle_res, funcIndex);
        double patch_relative_err = decayedRelativeError(patch_res, oracle_res, funcIndex);
        double ofp_relative_err = decayedRelativeError(ofp_res, oracle_res, funcIndex);

        metrics.decayed_origin_rel = max(metrics.decayed_origin_rel, origin_relative_err);
        metrics.decayed_aceso_rel = max(metrics.decayed_aceso_rel, patch_relative_err);
        metrics.decayed_ofp_repair_rel = max(metrics.decayed_ofp_repair_rel, ofp_relative_err);
        metrics.decayed_origin_abs = max(metrics.decayed_origin_abs, origin_absolute_err);
        metrics.decayed_aceso_abs = max(metrics.decayed_aceso_abs, patch_absolute_err);
        metrics.decayed_ofp_repair_abs = max(metrics.decayed_ofp_repair_abs, ofp_absolute_err);
    }

    return metrics;
}

RepairMetrics evaluateOFPMetrics(int funcIndex, double lbound, double rbound) {
    auto fixFunc = fixFuncList[funcIndex];
    auto oracleFunc = oracleFuncList[funcIndex];
    RepairMetrics metrics;

    auto inputs = randomOnValue(lbound, rbound, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = inputs[i];
        double ofp_res = fixFunc(xi);
        mpreal oracle_res = oracleFunc(xi);

        double ofp_absolute_err = decayedAbsoluteError(ofp_res, oracle_res);
        double ofp_relative_err = decayedRelativeError(ofp_res, oracle_res, funcIndex);

        metrics.stable_ofp_repair_rel = max(metrics.stable_ofp_repair_rel, ofp_relative_err);
        metrics.stable_ofp_repair_abs = max(metrics.stable_ofp_repair_abs, ofp_absolute_err);
    }

    double decayedHighMag = -2;
    if (hasPositiveLeftEndpointSingularity(funcIndex) || hasPositiveRightEndpointSingularity(funcIndex)) {
        decayedHighMag = log10(rbound - lbound);
    }
    auto mags = randomOnMagnitude(-20, decayedHighMag, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = (lbound + rbound)/2 + ((i%2==0) ? mags[i] : -mags[i]);
        if (hasPositiveLeftEndpointSingularity(funcIndex)) {
            xi = lbound + mags[i];
        } else if (hasPositiveRightEndpointSingularity(funcIndex)) {
            xi = rbound - mags[i];
        }
        double ofp_res = fixFunc(xi);
        mpreal oracle_res = oracleFunc(xi);

        double ofp_absolute_err = decayedAbsoluteError(ofp_res, oracle_res);
        double ofp_relative_err = decayedRelativeError(ofp_res, oracle_res, funcIndex);

        metrics.decayed_ofp_repair_rel = max(metrics.decayed_ofp_repair_rel, ofp_relative_err);
        metrics.decayed_ofp_repair_abs = max(metrics.decayed_ofp_repair_abs, ofp_absolute_err);
    }

    return metrics;
}

AutoRNPMetrics evaluateAutoRNPMetrics(int funcIndex, const std::function<double(double)>& autoRNPFunc,
                                      double lbound, double rbound) {
    auto originFunc = simpleFuncList[funcIndex];
    auto oracleFunc = oracleFuncList[funcIndex];
    AutoRNPMetrics metrics;

    auto inputs = randomOnValue(lbound, rbound, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = inputs[i];
        double origin_res = originFunc(xi);
        double autornp_res = autoRNPFunc(xi);
        mpreal oracle_res = oracleFunc(xi);

        double origin_absolute_err = decayedAbsoluteError(origin_res, oracle_res);
        double autornp_absolute_err = decayedAbsoluteError(autornp_res, oracle_res);
        double origin_relative_err = decayedRelativeError(origin_res, oracle_res, funcIndex);
        double autornp_relative_err = decayedRelativeError(autornp_res, oracle_res, funcIndex);

        metrics.stable_origin_abs = max(metrics.stable_origin_abs, origin_absolute_err);
        metrics.stable_autornp_abs = max(metrics.stable_autornp_abs, autornp_absolute_err);
        metrics.stable_origin_rel = max(metrics.stable_origin_rel, origin_relative_err);
        metrics.stable_autornp_rel = max(metrics.stable_autornp_rel, autornp_relative_err);
    }

    double decayedHighMag = -2;
    if (hasPositiveLeftEndpointSingularity(funcIndex) || hasPositiveRightEndpointSingularity(funcIndex)) {
        decayedHighMag = log10(rbound - lbound);
    }
    auto mags = randomOnMagnitude(-20, decayedHighMag, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = (lbound + rbound)/2 + ((i%2==0) ? mags[i] : -mags[i]);
        if (hasPositiveLeftEndpointSingularity(funcIndex)) {
            xi = lbound + mags[i];
        } else if (hasPositiveRightEndpointSingularity(funcIndex)) {
            xi = rbound - mags[i];
        }
        double origin_res = originFunc(xi);
        double autornp_res = autoRNPFunc(xi);
        mpreal oracle_res = oracleFunc(xi);

        double origin_absolute_err = decayedAbsoluteError(origin_res, oracle_res);
        double autornp_absolute_err = decayedAbsoluteError(autornp_res, oracle_res);
        double origin_relative_err = decayedRelativeError(origin_res, oracle_res, funcIndex);
        double autornp_relative_err = decayedRelativeError(autornp_res, oracle_res, funcIndex);

        metrics.decayed_origin_abs = max(metrics.decayed_origin_abs, origin_absolute_err);
        metrics.decayed_autornp_abs = max(metrics.decayed_autornp_abs, autornp_absolute_err);
        metrics.decayed_origin_rel = max(metrics.decayed_origin_rel, origin_relative_err);
        metrics.decayed_autornp_rel = max(metrics.decayed_autornp_rel, autornp_relative_err);
    }

    return metrics;
}

void writeAutoRNPHeader(std::ofstream& resultOut) {
    resultOut
        << "function_id,repeat_count,"
        << "stable_origin_abs,stable_autornp_abs,"
        << "stable_origin_rel,stable_autornp_rel,"
        << "decayed_origin_abs,decayed_autornp_abs,"
        << "decayed_origin_rel,decayed_autornp_rel\n";
}

void writeHerbieHeader(std::ofstream& resultOut) {
    resultOut
        << "function_id,repeat_count,"
        << "stable_origin_abs,stable_herbie_abs,"
        << "stable_origin_rel,stable_herbie_rel,"
        << "decayed_origin_abs,decayed_herbie_abs,"
        << "decayed_origin_rel,decayed_herbie_rel\n";
}

void writeMissingPatchRow(const char* resultPath, int funcIndex, int repeatCount,
                          void (*writeHeader)(std::ofstream&)) {
    bool writeHeaderNow = !fileExistsAndNotEmpty(resultPath);
    std::ofstream resultOut(resultPath, std::ios::app);
    if (writeHeaderNow) {
        writeHeader(resultOut);
    }
    resultOut << "\"" << functionIDs[funcIndex] << "\","
        << repeatCount
        << ",-,-,-,-,-,-,-,-" << std::endl;
}

int fitFullInfoAutoRNP(const AutoRNPEntry& entry, int repeatCount = 1) {
    const int funcIndex = entry.funcIndex;
    const char* resultPath = "result_autornp.csv";
    bool writeHeader = !fileExistsAndNotEmpty(resultPath);
    std::ofstream resultOut(resultPath, std::ios::app);
    if (writeHeader) {
        writeAutoRNPHeader(resultOut);
    }

    double lbound = errIntervals[funcIndex].first;
    double rbound = errIntervals[funcIndex].second;
    if (repeatCount < 1) {
        repeatCount = 1;
    }

    AutoRNPMetrics avgMetrics;
    for (int repeatIndex = 0; repeatIndex < repeatCount; repeatIndex++) {
        avgMetrics.add(evaluateAutoRNPMetrics(funcIndex, entry.patchFunc, lbound, rbound));
    }
    avgMetrics.divide((double)repeatCount);

    printf("Repeated Evaluation Count: %d\n", repeatCount);
    printf("Stable Area Improved By AutoRNP (abs avg): %.2e -> %.2e\n", avgMetrics.stable_origin_abs, avgMetrics.stable_autornp_abs);
    printf("Stable Area Improved By AutoRNP (rel avg): %.2e -> %.2e\n", avgMetrics.stable_origin_rel, avgMetrics.stable_autornp_rel);
    printf("Decayed Area Improved By AutoRNP (abs avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_abs, avgMetrics.decayed_autornp_abs);
    printf("Decayed Area Improved By AutoRNP (rel avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_rel, avgMetrics.decayed_autornp_rel);
    printf("\nEvaluation on Function: %s Done.\n", functionIDs[funcIndex].c_str());
    printf("==============\n");

    resultOut << "\"" << functionIDs[funcIndex] << "\","
        << repeatCount << ","
        << std::scientific
        << avgMetrics.stable_origin_abs << ","
        << avgMetrics.stable_autornp_abs << ","
        << avgMetrics.stable_origin_rel << ","
        << avgMetrics.stable_autornp_rel << ","
        << avgMetrics.decayed_origin_abs << ","
        << avgMetrics.decayed_autornp_abs << ","
        << avgMetrics.decayed_origin_rel << ","
        << avgMetrics.decayed_autornp_rel << std::endl;

    return 0;
}

int fitFullInfoHerbie(const HerbieEntry& entry, int repeatCount = 1) {
    const int funcIndex = entry.funcIndex;
    const char* resultPath = "result_herbie.csv";
    bool writeHeader = !fileExistsAndNotEmpty(resultPath);
    std::ofstream resultOut(resultPath, std::ios::app);
    if (writeHeader) {
        writeHerbieHeader(resultOut);
    }

    double lbound = errIntervals[funcIndex].first;
    double rbound = errIntervals[funcIndex].second;
    if (repeatCount < 1) {
        repeatCount = 1;
    }

    AutoRNPMetrics avgMetrics;
    for (int repeatIndex = 0; repeatIndex < repeatCount; repeatIndex++) {
        avgMetrics.add(evaluateAutoRNPMetrics(funcIndex, entry.patchFunc, lbound, rbound));
    }
    avgMetrics.divide((double)repeatCount);

    printf("Repeated Evaluation Count: %d\n", repeatCount);
    printf("Stable Area Improved By Herbie (abs avg): %.2e -> %.2e\n", avgMetrics.stable_origin_abs, avgMetrics.stable_autornp_abs);
    printf("Stable Area Improved By Herbie (rel avg): %.2e -> %.2e\n", avgMetrics.stable_origin_rel, avgMetrics.stable_autornp_rel);
    printf("Decayed Area Improved By Herbie (abs avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_abs, avgMetrics.decayed_autornp_abs);
    printf("Decayed Area Improved By Herbie (rel avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_rel, avgMetrics.decayed_autornp_rel);
    printf("\nEvaluation on Function: %s Done.\n", functionIDs[funcIndex].c_str());
    printf("==============\n");

    resultOut << "\"" << functionIDs[funcIndex] << "\","
        << repeatCount << ","
        << std::scientific
        << avgMetrics.stable_origin_abs << ","
        << avgMetrics.stable_autornp_abs << ","
        << avgMetrics.stable_origin_rel << ","
        << avgMetrics.stable_autornp_rel << ","
        << avgMetrics.decayed_origin_abs << ","
        << avgMetrics.decayed_autornp_abs << ","
        << avgMetrics.decayed_origin_rel << ","
        << avgMetrics.decayed_autornp_rel << std::endl;

    return 0;
}

void writeDifferentOrderHeader(std::ofstream& resultOut) {
    resultOut
        << "function_id,expansion_degree,repeat_count,"
        << "stable_ofp_repair_abs,stable_ofp_repair_rel,"
        << "decayed_ofp_repair_abs,decayed_ofp_repair_rel\n";
}

int fitDifferentOrderOFP(int funcIndex, int expansionDegree, int repeatCount = 10) {
    const char* resultPath = "result_differencet_order.csv";
    bool writeHeader = !fileExistsAndNotEmpty(resultPath);
    std::ofstream resultOut(resultPath, std::ios::app);
    if (writeHeader) {
        writeDifferentOrderHeader(resultOut);
    }

    double lbound = errIntervals[funcIndex].first;
    double rbound = errIntervals[funcIndex].second;
    if (repeatCount < 1) {
        repeatCount = 1;
    }

    set_expansion_terms(expansionDegree);
    RepairMetrics avgMetrics;
    for (int repeatIndex = 0; repeatIndex < repeatCount; repeatIndex++) {
        avgMetrics.add(evaluateOFPMetrics(funcIndex, lbound, rbound));
    }
    avgMetrics.divide((double)repeatCount);

    printf("OFP-Repair order %d, repeated evaluation count: %d\n", expansionDegree, repeatCount);
    printf("Stable Area OFP-Repair Error (abs avg): %.2e\n",
           avgMetrics.stable_ofp_repair_abs);
    printf("Stable Area OFP-Repair Error (rel avg): %.2e\n",
           avgMetrics.stable_ofp_repair_rel);
    printf("Decayed Area OFP-Repair Error (abs avg): %.2e\n",
           avgMetrics.decayed_ofp_repair_abs);
    printf("Decayed Area OFP-Repair Error (rel avg): %.2e\n",
           avgMetrics.decayed_ofp_repair_rel);
    printf("\nEvaluation on Function: %s Done.\n", functionIDs[funcIndex].c_str());
    printf("==============\n");

    resultOut << "\"" << functionIDs[funcIndex] << "\","
        << expansionDegree << ","
        << repeatCount << ","
        << std::scientific
        << avgMetrics.stable_ofp_repair_abs << ","
        << avgMetrics.stable_ofp_repair_rel << ","
        << avgMetrics.decayed_ofp_repair_abs << ","
        << avgMetrics.decayed_ofp_repair_rel << std::endl;

    return 0;
}

int fitFullInfo(int funcIndex, int weighted, int repeatCount = 1) {
    Serializer ser;
    ser.writeFuncIDToFile(functionIDs[funcIndex]);
    const char* resultPath = "result.csv";
    bool writeHeader = !fileExistsAndNotEmpty(resultPath);
    std::ofstream resultOut(resultPath, std::ios::app);
    if (writeHeader) {
        resultOut
            << "function_id,expansion_degree,repeat_count,"
            << "stable_origin_abs,stable_aceso_abs,stable_ofp_repair_abs,"
            << "stable_origin_rel,stable_aceso_rel,stable_ofp_repair_rel,"
            << "decayed_origin_abs,decayed_aceso_abs,decayed_ofp_repair_abs,"
            << "decayed_origin_rel,decayed_aceso_rel,decayed_ofp_repair_rel\n";
    }

    /////////////////// Synthesizing ///////////////////
    printf("Synthesizing Patch...\n");
    // measure time
    auto start_time = std::chrono::high_resolution_clock::now();

    Polynomial poly;
    double lbound = errIntervals[funcIndex].first;
    double rbound = errIntervals[funcIndex].second;
    int rootsNum  = defaultRootsNum;
    int highestDegree = defaultDegree;
    int polymode = defaultPoly;
    bool evaluateAceso = true;
    bool valid = false;

    evaluateAceso = !shouldSkipAcesoSynthesis(funcIndex);

    if (evaluateAceso) {
        Observation observation(funcIndex, lbound, rbound, rootsNum);
        observation.init();

        GSLFitModule gslfit;
        gslfit.init(highestDegree, weighted, polymode, lbound, rbound);
        gslfit.fitObservations(observation);
        poly = gslfit.getPoly();
        poly.checkacc0();

        auto end_time = std::chrono::high_resolution_clock::now();

        printf("Synthesized Polynomial: \n");
        poly.show();
        poly.writePolyToFile();
        int64_t time_on_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        printf("Time on Synthesizing: %ld ms\n\n", time_on_ms);

        /////////////////// Validation ///////////////////

        printf("Validating... ");
        start_time = std::chrono::high_resolution_clock::now();

        bool finitePoly = isFinitePolynomial(poly);
        int score = finitePoly ? observation.evaluatePoly(poly) : 0;
        if (finitePoly && score > 0) {
            printf("Valid.\n");
            valid = true;
        }
        else {
            printf("Invalid.\n");
            valid = false;
        }
        evaluateAceso = valid;
        ser.writeValidationToFile(valid);

        end_time = std::chrono::high_resolution_clock::now();

        time_on_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        printf("Time on Validation: %ld ms\n", time_on_ms);
        if (valid == false) {
            printf("[*DENIED Patch*] during validation on FuncIndex: %d\n", funcIndex);
            printf("[Notice]: Aceso patch metrics are set to inf for this function.\n");
        }
    }
    else {
        auto end_time = std::chrono::high_resolution_clock::now();
        int64_t time_on_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        printf("[Notice]: Aceso sampling is ignored for %s; patch metrics are set to inf.\n",
               functionIDs[funcIndex].c_str());
        printf("Time on Synthesizing: %ld ms\n", time_on_ms);
        ser.writeValidationToFile(false);
    }

    /////////////////// Evaluation with Oracle ///////////////////

    printf("\nEvaluation with Oracles...\n");
    if (repeatCount < 1) {
        repeatCount = 1;
    }
    RepairMetrics avgMetrics;
    for (int repeatIndex = 0; repeatIndex < repeatCount; repeatIndex++) {
        avgMetrics.add(evaluateRepairMetrics(funcIndex, poly, lbound, rbound, evaluateAceso));
    }
    avgMetrics.divide((double)repeatCount);

    printf("Repeated Evaluation Count: %d\n", repeatCount);
    printf("Stable Area Improved By Aceso (abs avg): %.2e -> %.2e\n", avgMetrics.stable_origin_abs, avgMetrics.stable_aceso_abs);
    printf("Stable Area Improved By Aceso (rel avg): %.2e -> %.2e\n", avgMetrics.stable_origin_rel, avgMetrics.stable_aceso_rel);
    printf("Stable Area Improved By OFP-Repair (abs avg): %.2e -> %.2e\n", avgMetrics.stable_origin_abs, avgMetrics.stable_ofp_repair_abs);
    printf("Stable Area Improved By OFP-Repair (rel avg): %.2e -> %.2e\n", avgMetrics.stable_origin_rel, avgMetrics.stable_ofp_repair_rel);
    printf("Decayed Area Improved By Aceso (abs avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_abs, avgMetrics.decayed_aceso_abs);
    printf("Decayed Area Improved By Aceso (rel avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_rel, avgMetrics.decayed_aceso_rel);
    printf("Decayed Area Improved By OFP-Repair (abs avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_abs, avgMetrics.decayed_ofp_repair_abs);
    printf("Decayed Area Improved By OFP-Repair (rel avg): %.2e -> %.2e\n", avgMetrics.decayed_origin_rel, avgMetrics.decayed_ofp_repair_rel);

    printf("\nEvaluation on Function: %s Done.\n", functionIDs[funcIndex].c_str());

    printf("==============\n");

    resultOut << "\"" << functionIDs[funcIndex] << "\","
        << get_expansion_terms() << ","
        << repeatCount << ","
        << std::scientific
        << avgMetrics.stable_origin_abs << ","
        << avgMetrics.stable_aceso_abs << ","
        << avgMetrics.stable_ofp_repair_abs << ","
        << avgMetrics.stable_origin_rel << ","
        << avgMetrics.stable_aceso_rel << ","
        << avgMetrics.stable_ofp_repair_rel << ","
        << avgMetrics.decayed_origin_abs << ","
        << avgMetrics.decayed_aceso_abs << ","
        << avgMetrics.decayed_ofp_repair_abs << ","
        << avgMetrics.decayed_origin_rel << ","
        << avgMetrics.decayed_aceso_rel << ","
        << avgMetrics.decayed_ofp_repair_rel << std::endl;

    ser.writeEvaluationToFile(
        avgMetrics.decayed_origin_rel,
        avgMetrics.decayed_aceso_rel,
        avgMetrics.decayed_origin_abs,
        avgMetrics.decayed_aceso_abs,
        avgMetrics.stable_origin_rel,
        avgMetrics.stable_aceso_rel,
        avgMetrics.stable_origin_abs,
        avgMetrics.stable_aceso_abs
    );
    return 0;
}

int fitFPResult(int funcIndex, int weighted) {
    auto originFunc = simpleFuncList[funcIndex];
    auto oracleFunc = oracleFuncList[funcIndex];

    Serializer ser;
    ser.writeFuncIDToFile(functionIDs[funcIndex]);
    // No Microstructure, cannot validate.
    ser.writeValidationToFile(true);

    /////////////////// Synthesizing ///////////////////
    Polynomial poly;
    double lbound = errIntervals[funcIndex].first;
    double rbound = errIntervals[funcIndex].second;
    int num  = defaultRootsNum;
    int highestDegree = defaultDegree;
    int polymode = defaultPoly;
    // Observation observation(funcIndex, lbound, rbound, num);
    // observation.init();

    GSLFitModule gslfit;
    // Fit with direct FP values
    gslfit.init(highestDegree, weighted, polymode, lbound, rbound);
    std::vector<double> xs;
    for (int k = 0; k <= num; k++) {
        double fpk = k;
        double fpn = num+1;
        double rt = cosl((2*fpk+1)/(2*fpn)*myPi);
        rt = (lbound+rbound)/2 + rt*(rbound-lbound)/2;
        xs.push_back(rt);
    }
    std::vector<double> ys;
    std::vector<double> rxs;
    poly.lbound = lbound;
    poly.rbound = rbound;
    for (double xi : xs) {
        double y = originFunc(xi);
        if (isnan(y) || isinf(y))
            continue;
        else {
            ys.push_back(originFunc(xi));
            rxs.push_back(poly.reduction(xi));
        }
    }

    std::vector<double> coef = gslfit._fitCenterPointsOnly(rxs, ys);
    poly.coef = coef;
    poly.mode = defaultPoly;

    printf("Synthesized Polynomial: \n");
    poly.show();

    /////////////////// Evaluation with Oracle ///////////////////
    printf("\nEvaluation with Oracles...\n");

    // On accurate area
    double stable_origin_max_relative_err = 0;
    double stable_patch_max_relative_err = 0;

    double stable_origin_max_absolute_err = 0;
    double stable_patch_max_absolute_err = 0;

    auto inputs = randomOnValue(lbound, rbound, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = inputs[i];
        double origin_res = originFunc(xi);
        double patch_res = poly.evaluate(xi);
        mpreal oracle_res = oracleFunc(xi);

        // printf("=============\nx: %.3e\n", xi);
        // printf("origin: %.15e\n", origin_res);
        // printf("poly:   %.15e\n", patch_res);
        // printf("oracle: %.15e\n", (double)(oracle_res));;
        // getchar();

        double origin_absolute_err = (double)fabs(origin_res - oracle_res);
        double patch_absolute_err  = (double)fabs(patch_res - oracle_res);
        double origin_relative_err = (double)fabs((origin_res - oracle_res) / oracle_res);
        double patch_relative_err  = (double)fabs((patch_res - oracle_res) / oracle_res);

        stable_origin_max_relative_err = max(stable_origin_max_relative_err, origin_relative_err);
        stable_patch_max_relative_err  = max(stable_patch_max_relative_err, patch_relative_err);
        stable_origin_max_absolute_err = max(stable_origin_max_absolute_err, origin_absolute_err);
        stable_patch_max_absolute_err  = max(stable_patch_max_absolute_err, patch_absolute_err);
    }
    printf("Stable Area Improved By Naive Fitting (abs): %.2e -> %.2e\n", stable_origin_max_absolute_err, stable_patch_max_absolute_err);
    printf("Stable Area Improved By Naive Fitting (rel): %.2e -> %.2e\n", stable_origin_max_relative_err, stable_patch_max_relative_err);

    // On decayed area
    double decayed_origin_max_relative_err = 0;
    double decayed_patch_max_relative_err = 0;

    double decayed_origin_max_absolute_err = 0;
    double decayed_patch_max_absolute_err = 0;

    auto mags = randomOnMagnitude(-20, -2, evaluationNum);
    for (int i = 0; i < evaluationNum; i++) {
        double xi = (lbound + rbound)/2 + ((i%2==0) ? mags[i] : -mags[i]);
        double origin_res = originFunc(xi);
        double patch_res = poly.evaluate(xi);
        mpreal oracle_res = oracleFunc(xi);

        double origin_absolute_err = (double)fabs(origin_res - oracle_res);
        double patch_absolute_err  = (double)fabs(patch_res - oracle_res);
        double origin_relative_err = (double)fabs((origin_res - oracle_res) / oracle_res);
        double patch_relative_err  = (double)fabs((patch_res - oracle_res) / oracle_res);

        decayed_origin_max_relative_err = max(decayed_origin_max_relative_err, origin_relative_err);
        decayed_patch_max_relative_err  = max(decayed_patch_max_relative_err, patch_relative_err);
        decayed_origin_max_absolute_err = max(decayed_origin_max_absolute_err, origin_absolute_err);
        decayed_patch_max_absolute_err  = max(decayed_patch_max_absolute_err, patch_absolute_err);
    }
    printf("Decayed Area Improved By Naive Fitting (abs): %.2e -> %.2e\n", decayed_origin_max_absolute_err, decayed_patch_max_absolute_err);
    printf("Decayed Area Improved By Naive Fitting (rel): %.2e -> %.2e\n", decayed_origin_max_relative_err, decayed_patch_max_relative_err);

    printf("\nEvaluation on Function: %s Done.\n", functionIDs[funcIndex].c_str());

    printf("==============\n");

    ser.writeEvaluationToFile(
        decayed_origin_max_relative_err,
        decayed_patch_max_relative_err,
        decayed_origin_max_absolute_err,
        decayed_patch_max_absolute_err,
        stable_origin_max_relative_err,
        stable_patch_max_relative_err,
        stable_origin_max_absolute_err,
        stable_patch_max_absolute_err
    );
    return 0;
}

void testFunc(int funcIndex) {
    auto originFunc = simpleFuncList[funcIndex];
    auto oracleFunc = oracleFuncList[funcIndex];

    printf("Testing Function Number: %d\n", funcIndex);
    while (1) {
        double x;
        printf("Input >> ");
        scanf("%lf", &x);
        printf("origin: %.15e\n", originFunc(x));
        printf("oracle: %.15e\n", (double)oracleFunc(x));
    }
    return;
}

int matchFuncID(char* ID) {
    int singleindex=-1;
    for (int i = 0; i < simpleFuncList.size(); i++) {
        if (strcmp(ID, functionIDs[i].c_str()) == 0) {
            return i;
        }
    }
    return -1;
}

bool isIntegerArg(const char* s) {
    if (s == nullptr || *s == '\0') {
        return false;
    }
    if (*s == '-' || *s == '+') {
        ++s;
    }
    if (*s == '\0') {
        return false;
    }
    while (*s != '\0') {
        if (*s < '0' || *s > '9') {
            return false;
        }
        ++s;
    }
    return true;
}

bool fileExistsAndNotEmpty(const char* path) {
    struct stat st;
    return stat(path, &st) == 0 && st.st_size > 0;
}

int main(int argc, char** argv) {

    Serializer ser;
    ser.cleanFuncIDFile();
    ser.cleanValidationFile();
    ser.cleanEvaluationFile();
    ser.cleanInfoLevelFile();

    if (argc > 1) {
        if (strcmp(argv[1],"fullinfo_ofp")==0) {
            ser.writeInfoLevelToFile(FullInfoLevel);
            int repeatCount = 1;
            if (argc == 2 || (argc >= 3 && isIntegerArg(argv[2]))) {
                if (argc >= 3) {
                    set_expansion_terms(atoi(argv[2]));
                }
                if (argc >= 4) {
                    repeatCount = atoi(argv[3]);
                }
                for (int i = 0; i < simpleFuncList.size(); i++) {
                    fitFullInfo(i, USE_WEIGHT, repeatCount);
                }
            }
            else if (argc >= 3) {
                repeatCount = 1;
                if (argc >= 4) {
                    set_expansion_terms(atoi(argv[3]));
                }
                if (argc >= 5) {
                    repeatCount = atoi(argv[4]);
                }
                int funcIndex = matchFuncID(argv[2]);
                if (funcIndex == -1) {
                    printf("bin/testground.out fullinfo_ofp [ expansionDegree ] [ repeatCount ]\n");
                    printf("bin/testground.out fullinfo_ofp [ functionID ] [ expansionDegree ] [ repeatCount ]\n");
                    return 0;
                }
                fitFullInfo(funcIndex, USE_WEIGHT, repeatCount);
            }
        }
        else if (strcmp(argv[1],"fullinfo_autornp")==0) {
            int repeatCount = 1;
            if (argc >= 3) {
                repeatCount = atoi(argv[2]);
            }
            for (int i = 0; i < simpleFuncList.size(); i++) {
                const AutoRNPEntry* entry = findAutoRNPEntry(i);
                if (entry == nullptr) {
                    writeMissingPatchRow("result_autornp.csv", i, repeatCount, writeAutoRNPHeader);
                } else {
                    fitFullInfoAutoRNP(*entry, repeatCount);
                }
            }
        }
        else if (strcmp(argv[1],"fullinfo_herbie")==0) {
            int repeatCount = 1;
            if (argc >= 3) {
                repeatCount = atoi(argv[2]);
            }
            for (int i = 0; i < simpleFuncList.size(); i++) {
                const HerbieEntry* entry = findHerbieEntry(i);
                if (entry == nullptr) {
                    writeMissingPatchRow("result_herbie.csv", i, repeatCount, writeHerbieHeader);
                } else {
                    fitFullInfoHerbie(*entry, repeatCount);
                }
            }
        }
        else if (strcmp(argv[1],"fullinfo_differencet_order")==0) {
            int maxExpansionDegree = 10;
            int repeatCount = 10;
            if (argc >= 3) {
                maxExpansionDegree = atoi(argv[2]);
            }
            if (argc >= 4) {
                repeatCount = atoi(argv[3]);
            }
            if (maxExpansionDegree < 1) {
                maxExpansionDegree = 1;
            }
            for (int degree = 1; degree <= maxExpansionDegree; degree++) {
                for (int i = 0; i < simpleFuncList.size(); i++) {
                    fitDifferentOrderOFP(i, degree, repeatCount);
                }
            }
        }
        else
        if (strcmp(argv[1],"noinfo")==0) {
            if (argc == 2) {
                ser.writeInfoLevelToFile(NoInfoLevel);
                for (int i = 0; i < simpleFuncList.size(); i++) {
                    // Useful in evaluation part.
                    fitFPResult(i, NO_WEIGHT); // Fitting with no information
                }
            }
            else {
                int funcIndex = matchFuncID(argv[2]);
                if (funcIndex == -1) {
                    printf("bin/testground.out noinfo [ functionID ]\n");
                    return 0;
                }
                else {
                    fitFPResult(funcIndex, NO_WEIGHT);
                }
            }
        }
        else if (strcmp(argv[1],"fullinfo")==0) {
            if (argc == 2) {
                ser.writeInfoLevelToFile(FullInfoLevel);
                for (int i = 0; i < simpleFuncList.size(); i++) {
                    // Useful in evaluation part.
                    fitFullInfo(i, USE_WEIGHT); // Fitting with full information from micro-structure
                }
            }
            else {
                int funcIndex = matchFuncID(argv[2]);
                if (funcIndex == -1) {
                    printf("bin/testground.out fullinfo [ functionID ]\n");
                    return 0;
                }
                else {
                    fitFullInfo(funcIndex, USE_WEIGHT);
                }
            }
        }
        else if (strcmp(argv[1],"onlyweight")==0) {
            ser.writeInfoLevelToFile(1);
            for (int i = 0; i < simpleFuncList.size(); i++) {
                // Useful in evaluation part.
                fitFPResult(i, USE_WEIGHT); // Fitting with full information from micro-structure
            }
        }
        else if (strcmp(argv[1],"runtime")==0) {
            for (int i = 0; i < simpleFuncList.size(); i++) {
                measurePatchRuntime(i, USE_WEIGHT);
            }
        }
        else {
            printf("bin/testground.out [ fullinfo / noinfo ]\n");
        }
        return 0;
    }

    // Default: Full info
    ser.writeInfoLevelToFile(FullInfoLevel);
    for (int i = 0; i < simpleFuncList.size(); i++) {
        // Useful in evaluation part.
        fitFullInfo(i, USE_WEIGHT); // Fitting with full information from micro-structure
    }

    // for (int i = 0; i < simpleFuncList.size(); i++) {
    //     // Unnecessary to provide figure, just words.
    //     // fitFPResult(i, USE_WEIGHT); // Fitting with FP, with weight
    //     // fitFullInfo(i, NO_WEIGHT);  // Fitting with mean, without weight
    //     // getchar();
    // }
    // fitFullInfo(7, USE_WEIGHT);

    return 0;
}
