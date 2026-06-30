#pragma once

#include <gsl/gsl_multifit.h>
#include <vector>
#include <cassert>

std::pair<double, double> fitPower2(const std::vector<double>& xs, const std::vector<double>& ys) {
	assert(xs.size() == ys.size());
	int rootsNum = xs.size();
	// only need for c1 and c2;
	int paramsNum = 2;
	double xi, yi, ei;
	double xij;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;
	double chisq;
	// Allocate
	// X, usually xi, xi^2, xi^3... (Expansion polynomial)
	X = gsl_matrix_alloc(rootsNum, paramsNum);
	// y, the observation of y value
	y = gsl_vector_alloc(rootsNum);
	c = gsl_vector_alloc(paramsNum);
	cov = gsl_matrix_alloc(paramsNum, paramsNum);
	// Feed data into X, y
	for (int i = 0; i < rootsNum; i++) {
		// Set X
		xi = xs[i];
		gsl_matrix_set(X, i, 0, xi);
		gsl_matrix_set(X, i, 1, xi*xi);
		// Set y
		yi = ys[i];
		gsl_vector_set(y, i, yi);
	}
	// Fitting using GSL's method
	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(rootsNum, paramsNum);
	int status = gsl_multifit_linear(X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free(work);
	double c1, c2;
	if (status == 0) {
		c1 = gsl_vector_get(c, 0);
		c2 = gsl_vector_get(c, 1);
	}
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

	return std::make_pair(c1, c2);
}
