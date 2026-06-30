#ifndef TARGETS_H
#define TARGETS_H

#include <functional>
#include <vector>

extern "C" void set_expansion_terms(int terms);
extern "C" int get_expansion_terms();

// targets.c
extern "C" double exp_BI(double);
extern "C" double bJ_sin(double);
extern "C" double di_tan(double);
extern "C" double log_erf(double);
extern "C" double acos_fd(double);
extern "C" double ei(double);
extern "C" double Q1_W(double);
extern "C" double bj_tan(double);
extern "C" double Si_tan(double);
extern "C" double by_psi(double);
extern "C" double fdm_log(double);
extern "C" double eQ_sqrt(double);
extern "C" double W_var(double);
extern "C" double W_log(double);
extern "C" double cos_x2(double);
extern "C" double exp_2(double);
extern "C" double cos_sin(double);
extern "C" double sin_sin(double);
extern "C" double tan_tan(double);
extern "C" double cos_cos(double);
extern "C" double exp_exp(double);
extern "C" double exp_1(double);
extern "C" double x_tan(double);
extern "C" double log_log(double);
extern "C" double log_x(double);
extern "C" double sqrt_exp(double);
extern "C" double sin_tan(double);
extern "C" double exp_x(double);
extern "C" double x_x2(double);
extern "C" double sci3545_2(double);
extern "C" double sci4034_1(double);
extern "C" double sci4034_2(double);
extern "C" double sci3547(double);
extern "C" double sci3545_1(double);
extern "C" double pow_df(double);
extern "C" double chi_ci(double);
extern "C" double fc_bj(double);
extern "C" double gb_sqrt(double);
extern "C" double l1_l2(double);
extern "C" double hyp_g2(double);
extern "C" double sympy_28514_log_ratio(double);
extern "C" double sympy_12671_log_cosh(double);
extern "C" double statsmodels_1604_logistic(double);
extern "C" double statsmodels_826_cloglog(double);
extern "C" double scipy_17194_expm1_kernel(double);
extern "C" double pytorch_39242_log1mexp(double);


extern "C" double exp_BI_ofp_fix(double);
extern "C" double bJ_sin_ofp_fix(double);
extern "C" double di_tan_ofp_fix(double);
extern "C" double log_erf_ofp_fix(double);
extern "C" double acos_fd_ofp_fix(double);
extern "C" double ei_ofp_fix(double);
extern "C" double Q1_W_ofp_fix(double);
extern "C" double bj_tan_ofp_fix(double);
extern "C" double Si_tan_ofp_fix(double);
extern "C" double by_psi_ofp_fix(double);
extern "C" double fdm_log_ofp_fix(double);
extern "C" double eQ_sqrt_ofp_fix(double);
extern "C" double W_var_ofp_fix(double);
extern "C" double W_log_ofp_fix(double);
extern "C" double cos_x2_ofp_fix(double);
extern "C" double exp_2_ofp_fix(double);
extern "C" double cos_sin_ofp_fix(double);
extern "C" double sin_sin_ofp_fix(double);
extern "C" double tan_tan_ofp_fix(double);
extern "C" double cos_cos_ofp_fix(double);
extern "C" double exp_exp_ofp_fix(double);
extern "C" double exp_1_ofp_fix(double);
extern "C" double x_tan_ofp_fix(double);
extern "C" double log_log_ofp_fix(double);
extern "C" double log_x_ofp_fix(double);
extern "C" double sqrt_exp_ofp_fix(double);
extern "C" double sin_tan_ofp_fix(double);
extern "C" double exp_x_ofp_fix(double);
extern "C" double x_x2_ofp_fix(double);
extern "C" double sci3545_2_ofp_fix(double);
extern "C" double sci4034_1_ofp_fix(double);
extern "C" double sci4034_2_ofp_fix(double);
extern "C" double sci3547_ofp_fix(double);
extern "C" double sci3545_1_ofp_fix(double);
extern "C" double pow_df_ofp_fix(double);
extern "C" double chi_ci_ofp_fix(double);
extern "C" double fc_bj_ofp_fix(double);
extern "C" double gb_sqrt_ofp_fix(double);
extern "C" double l1_l2_ofp_fix(double);
extern "C" double hyp_g2_ofp_fix(double);
extern "C" double sympy_28514_log_ratio_ofp_fix(double);
extern "C" double sympy_12671_log_cosh_ofp_fix(double);
extern "C" double statsmodels_1604_logistic_ofp_fix(double);
extern "C" double statsmodels_826_cloglog_ofp_fix(double);
extern "C" double scipy_17194_expm1_kernel_ofp_fix(double);
extern "C" double pytorch_39242_log1mexp_ofp_fix(double);

const std::vector<std::function<double(double)>> simpleFuncList = {
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
    sci3545_2,
    sci4034_1,
    sci4034_2,
    sci3547,
    sci3545_1,
    sympy_28514_log_ratio,
    sympy_12671_log_cosh,
    statsmodels_1604_logistic,
    statsmodels_826_cloglog,
    scipy_17194_expm1_kernel,
    pytorch_39242_log1mexp,
};

const std::vector<std::function<double(double)>> fixFuncList = {
    exp_BI_ofp_fix,
    bJ_sin_ofp_fix,
    di_tan_ofp_fix,
    log_erf_ofp_fix,
    acos_fd_ofp_fix,
    ei_ofp_fix,
    Q1_W_ofp_fix,
    bj_tan_ofp_fix,
    Si_tan_ofp_fix,
    by_psi_ofp_fix,
    fdm_log_ofp_fix,
    eQ_sqrt_ofp_fix,
    W_var_ofp_fix,
    W_log_ofp_fix,
    pow_df_ofp_fix,
    chi_ci_ofp_fix,
    fc_bj_ofp_fix,
    cos_x2_ofp_fix,
    exp_2_ofp_fix,
    cos_sin_ofp_fix,
    sin_sin_ofp_fix,
    tan_tan_ofp_fix,
    cos_cos_ofp_fix,
    exp_exp_ofp_fix,
    exp_1_ofp_fix,
    x_tan_ofp_fix,
    log_log_ofp_fix,
    log_x_ofp_fix,
    sqrt_exp_ofp_fix,
    sin_tan_ofp_fix,
    exp_x_ofp_fix,
    x_x2_ofp_fix,
    sci3545_2_ofp_fix,
    sci4034_1_ofp_fix,
    sci4034_2_ofp_fix,
    sci3547_ofp_fix,
    sci3545_1_ofp_fix,
    sympy_28514_log_ratio_ofp_fix,
    sympy_12671_log_cosh_ofp_fix,
    statsmodels_1604_logistic_ofp_fix,
    statsmodels_826_cloglog_ofp_fix,
    scipy_17194_expm1_kernel_ofp_fix,
    pytorch_39242_log1mexp_ofp_fix,
};

#endif
