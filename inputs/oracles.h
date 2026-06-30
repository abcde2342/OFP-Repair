#pragma once

#include <mpreal.h>

using mpfr::mpreal;

mpreal oracle_exp_BI(double x);
mpreal oracle_bJ_sin(double x);
mpreal oracle_di_tan(double x);
mpreal oracle_log_erf(double x);
mpreal oracle_acos_fd(double x);
mpreal oracle_ei(double x);
mpreal oracle_Q1_W(double x);
mpreal oracle_bj_tan(double x);
mpreal oracle_Si_tan(double x);
mpreal oracle_by_psi(double x);
mpreal oracle_fdm_log(double x);
mpreal oracle_eQ_sqrt(double x);
mpreal oracle_W_var(double x);
mpreal oracle_W_log(double x);
mpreal oracle_cos_x2(double x);
mpreal oracle_exp_2(double x);
mpreal oracle_cos_sin(double x);
mpreal oracle_sin_sin(double x);
mpreal oracle_tan_tan(double x);
mpreal oracle_cos_cos(double x);
mpreal oracle_exp_exp(double x);
mpreal oracle_exp_1(double x);
mpreal oracle_x_tan(double x);
mpreal oracle_log_log(double x);
mpreal oracle_log_x(double x);
mpreal oracle_sqrt_exp(double x);
mpreal oracle_sin_tan(double x);
mpreal oracle_exp_x(double x);
mpreal oracle_x_x2(double x);
mpreal oracle_sci3545_2(double x);
mpreal oracle_sci4034_1(double x);
mpreal oracle_sci4034_2(double x);
mpreal oracle_sci3547(double x);
mpreal oracle_sci3545_1(double x);
mpreal oracle_pow_df(double x);
mpreal oracle_chi_ci(double x);
mpreal oracle_fc_bj(double x);
mpreal oracle_sympy_28514_log_ratio(double x);
mpreal oracle_sympy_12671_log_cosh(double x);
mpreal oracle_statsmodels_1604_logistic(double x);
mpreal oracle_statsmodels_826_cloglog(double x);
mpreal oracle_scipy_17194_expm1_kernel(double x);
mpreal oracle_pytorch_39242_log1mexp(double x);

const std::vector<std::function<mpreal(double)>> oracleFuncList = {
    oracle_exp_BI,
    oracle_bJ_sin,
    oracle_di_tan,
    oracle_log_erf,
    oracle_acos_fd,
    oracle_ei,
    oracle_Q1_W,
    oracle_bj_tan,
    oracle_Si_tan,
    oracle_by_psi,
    oracle_fdm_log,
    oracle_eQ_sqrt,
    oracle_W_var,
    oracle_W_log,
    oracle_pow_df,
    oracle_chi_ci,
    oracle_fc_bj,
    oracle_cos_x2,
    oracle_exp_2,
    oracle_cos_sin,
    oracle_sin_sin,
    oracle_tan_tan,
    oracle_cos_cos,
    oracle_exp_exp,
    oracle_exp_1,
    oracle_x_tan,
    oracle_log_log,
    oracle_log_x,
    oracle_sqrt_exp,
    oracle_sin_tan,
    oracle_exp_x,
    oracle_x_x2,
    oracle_sci3545_2,
    oracle_sci4034_1,
    oracle_sci4034_2,
    oracle_sci3547,
    oracle_sci3545_1,
    oracle_sympy_28514_log_ratio,
    oracle_sympy_12671_log_cosh,
    oracle_statsmodels_1604_logistic,
    oracle_statsmodels_826_cloglog,
    oracle_scipy_17194_expm1_kernel,
    oracle_pytorch_39242_log1mexp,
};
