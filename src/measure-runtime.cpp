#include "../patches/patches.h"

#include <math.h>
#include <stdio.h>

#include <chrono>
#include <random>
#include <vector>

typedef double (*undertest)(double);

double test_runtime(const char* id, int nums, double lbound, double rbound, undertest patched, undertest origin, bool valid) {
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(lbound, rbound);
  double d;

  std::chrono::high_resolution_clock::time_point t1, t2;
  std::chrono::nanoseconds ns1, ns2;
  double avg_ns1, avg_ns2;

  std::vector<double> v(nums);
  for (int i = 0; i < nums; i++) {
    v[i] = dist(e2);
  }

  double res1 = 0.0;
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < nums; i++) {
    d = v[i];
    res1 += patched(d);
  }
  t2 = std::chrono::high_resolution_clock::now();
  ns1 = std::chrono::duration_cast<std::chrono::nanoseconds> (t2-t1);
  avg_ns1 = (double) ns1.count() / nums;

  double res2 = 0.0;
  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < nums; i++) {
    d = v[i];
    res2 += origin(d);
  }
  t2 = std::chrono::high_resolution_clock::now();
  ns2 = std::chrono::duration_cast<std::chrono::nanoseconds> (t2-t1);
  avg_ns2 = (double) ns2.count() / nums;

  if (valid)
    printf("ID: %s\t\t%.2f\t\t%.2f\n", id, avg_ns2, avg_ns1);
  else
    printf("ID: %s\t\t%.2f\t\t%s\n", id, avg_ns2, "n/a");
  // printf("%.5e\n", res1);
  // printf("%.5e\n", res2);
  return res1 + res2;
}

int main() {
  printf("Time for origin / patched (nanoseconds)\n");
  test_runtime("1 exp_BI", 1000000, -0.01, 0.01, patched_exp_BI, exp_BI, 1);
  test_runtime("2 bJ_sin", 1000000, -0.01, 0.01, patched_bJ_sin, bJ_sin, 1);
  test_runtime("3 di_tan", 1000000, -0.01, 0.01, patched_di_tan, di_tan, 1);
  test_runtime("4 log_erf", 1000000, -0.01, 0.01, patched_log_erf, log_erf, 1);
  test_runtime("5 acos_fd", 1000000, -0.01, 0.01, patched_acos_fd, acos_fd, 1);
  test_runtime("6 ei", 1000000, -0.01, 0.01, patched_ei, ei, 1);
  test_runtime("7 Q1_W", 1000000, -0.01, 0.01, patched_Q1_W, Q1_W, 1);
  test_runtime("8 bj_tan", 1000000, -0.01, 0.01, patched_bj_tan, bj_tan, 1);
  test_runtime("9 Si_tan", 1000000, -0.01, 0.01, patched_Si_tan, Si_tan, 1);
  test_runtime("10 by_psi", 1000000, -0.01, 0.01, patched_by_psi, by_psi, 1);
  test_runtime("11 fdm_log", 1000000, -0.01, 0.01, patched_fdm_log, fdm_log, 1);
  test_runtime("12 eQ_sqrt", 1000000, -0.01, 0.01, patched_eQ_sqrt, eQ_sqrt, 1);
  test_runtime("13 W_var", 1000000, 2.7082818284590453, 2.728281828459045, patched_W_var, W_var, 1);
  test_runtime("14 W_log", 1000000, 2.7082818284590453, 2.728281828459045, patched_W_log, W_log, 1);
  test_runtime("15 pow_df", 1000000, -0.01, 0.01, patched_pow_df, pow_df, 1);
  test_runtime("16 chi_ci", 1000000, -0.01, 0.01, patched_chi_ci, chi_ci, 1);
  test_runtime("17 fc_bj", 1000000, -0.01, 0.01, patched_fc_bj, fc_bj, 1);
  test_runtime("18 gb_sqrt", 1000000, 2.99, 3.01, patched_gb_sqrt, gb_sqrt, 1);
  test_runtime("19 l1_l2", 1000000, 0.99, 1.01, patched_l1_l2, l1_l2, 1);
  test_runtime("20 hyp_g2", 1000000, 0.99, 1.01, patched_hyp_g2, hyp_g2, 1);
  test_runtime("s1 cos_x2", 1000000, -0.01, 0.01, patched_cos_x2, cos_x2, 1);
  test_runtime("s2 exp_2", 1000000, -0.01, 0.01, patched_exp_2, exp_2, 1);
  test_runtime("s3 cos_sin", 1000000, -0.01, 0.01, patched_cos_sin, cos_sin, 1);
  test_runtime("s4 sin_sin", 1000000, 1.5607963267948965, 1.5807963267948966, patched_sin_sin, sin_sin, 1);
  test_runtime("s5 tan_tan", 1000000, -0.01, 0.01, patched_tan_tan, tan_tan, 1);
  test_runtime("s6 cos_cos", 1000000, -0.01, 0.01, patched_cos_cos, cos_cos, 1);
  test_runtime("s7 exp_exp", 1000000, -0.01, 0.01, patched_exp_exp, exp_exp, 1);
  test_runtime("s8 exp_1", 1000000, -0.01, 0.01, patched_exp_1, exp_1, 1);
  test_runtime("s9 x_tan", 1000000, -0.01, 0.01, patched_x_tan, x_tan, 1);
  test_runtime("s10 log_log", 1000000, -0.01, 0.01, patched_log_log, log_log, 1);
  test_runtime("s11 log_x", 1000000, -0.01, 0.01, patched_log_x, log_x, 1);
  test_runtime("s12 sqrt_exp", 1000000, -0.01, 0.01, patched_sqrt_exp, sqrt_exp, 1);
  test_runtime("s13 sin_tan", 1000000, -0.01, 0.01, patched_sin_tan, sin_tan, 1);
  test_runtime("s14 exp_x", 1000000, -0.01, 0.01, patched_exp_x, exp_x, 1);
  test_runtime("s15 x_x2", 1000000, 0.99, 1.01, patched_x_x2, x_x2, 1);
  test_runtime("c1 sci3545_2", 1000000, 1e-6, 0.01, patched_control_1, sci3545_2, 0);
  test_runtime("c2 control_2", 1000000, -0.01, 0.01, patched_control_2, control_2, 0);
  test_runtime("c3 control_3", 1000000, -0.01, 0.01, patched_control_3, control_3, 0);
  test_runtime("c4 sci3547", 1000000, -0.01, 0.01, patched_control_4, sci3547, 0);
}
