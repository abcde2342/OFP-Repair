#ifndef GSLFIT_H
#define GSLFIT_H

#include <gsl/gsl_multifit.h>
#include <cstdio>
#include <vector>

#include "observation.h"
#include "serializer.h"
#include "local-observation.h"
#include "polynomial.h"

class GSLFitModule {
private:
	// Observations
	std::vector<LocalObservation> rootsInfo;

	// highestDegree
	int highestDegree = 6;
	// with weight or not
	int weighted = 1;

	// Results part
	// Params of polynomial
	int status = -1;
	
	Polynomial poly;
	std::vector<std::vector<double> > xCov;
	double chisq;

private:
	// double expansion_based_xij(double xi, int j) {
	// 	// return pow(xi, j)
	// 	if (j < 0) return expansion_based_xij(xi, -j);
	// 	if (j == 0) return 1.0;
	// 	if (j == 1) return xi;

	// 	if ((j & 1) == 0) {
	// 		double half = expansion_based_xij(xi, j/2);
	// 		return half*half;
	// 	}
	// 	else {
	// 		double half = expansion_based_xij(xi, j/2);
	// 		return xi * half * half;
	// 	}
	// }

public:
	void init(int degree, int wtd, int polymode, double lb, double rb) {
		highestDegree = degree;
		weighted = wtd;
		poly.mode = polymode;
		poly.lbound = lb;
		poly.rbound = rb;

		// Force polyBaseMode = 0;
		if (poly.mode != 0 && poly.mode != 1) {
			printf("[src/git-fit.h] : this mode is not implement yet!\n");
			poly.mode = 0;
		}
	}

	Polynomial getPoly() {
		if (status != 0) {
			printf("Fitting Result not available.\n");
			return Polynomial(0,0,0);
		}
		return poly;
	}

	std::vector< std::vector<double> > getCov() {
		if (status != 0) {
			printf("Cov Result not available.\n");
			return std::vector< std::vector<double> >{};
		}
		return xCov;
	}

	double getChisq() {
		if (status != 0) {
			printf("chisq not not available.\n");
			return 0;
		}
		return chisq;
	}

	int fitObservations(const Observation &observation) {
		/*
		* Ref: https://www.gnu.org/software/gsl/doc/html/lls.html#multi-parameter-linear-regression-example
		*
		* By define different X, we can use Expansion or Chebyshev polynomials
		*/

		// Prepare Parameters
		int rootsNum = observation.getSize();
		int paramsNum = highestDegree + 1;

	 	// Variables
		double xi, yi, ei, vari;
		double xij;
		gsl_matrix *X, *cov;
		gsl_vector *y, *w, *c;

		// Allocate
		// X, usually xi, xi^2, xi^3... (Expansion polynomial)
		X = gsl_matrix_alloc(rootsNum, paramsNum);
		// y, the observation of y value
		y = gsl_vector_alloc(rootsNum);
		// weight, usually the 1 / (err^2)
		w = gsl_vector_alloc(rootsNum);

		// store the results of polynomial parameters
		c = gsl_vector_alloc(paramsNum);
		cov = gsl_matrix_alloc(paramsNum, paramsNum);

		// Feed data into X, y, and w
		for (int i = 0; i < rootsNum; i++) {
			LocalObservation rootInfo = observation.getObservation(i);

			// Set X Vandermonde matrix
			xi = rootInfo.xMean;
			auto vanderVector = poly.vandermonde(xi, paramsNum);
			for (int j = 0; j < paramsNum; j++) {
				xij = vanderVector[j];
				// xij = poly.calcOnTerm(xi, j);
				gsl_matrix_set(X, i, j, xij);
			}

			// Set y
			yi = rootInfo.yMean;
			gsl_vector_set(y, i, yi);

			// Set weight of error
			// ei = (rootInfo.yMax - rootInfo.yMin) / 2.0;
			vari = rootInfo.yVar;
			if (weighted) {
				// gsl_vector_set(w, i, 1.0/(ei*ei));
				gsl_vector_set(w, i, 1/(1e-30+vari));
			}
			else {
				gsl_vector_set(w, i, 1.0);
			}
		}

		// Fitting using GSL's method
		gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(rootsNum, paramsNum);
		status = gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
		gsl_multifit_linear_free(work);

		if (status == 0) {

			// Store params from C to C++
			poly.coef.resize(paramsNum);
			for (int i = 0; i < paramsNum; i++) {
				poly.coef[i] = gsl_vector_get(c, i);
			}

			// Store params from C to C++
			xCov.resize(paramsNum);
			for (int i = 0; i < paramsNum; i++) {
				xCov[i].resize(paramsNum);
				for (int j = 0; j < paramsNum; j++) {
					xCov[i][j] = gsl_matrix_get(cov, i, j);
				}
			}
		}

		// Free space
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);

		return status;
	}

	std::vector<double> _fitCenterPointsOnly(std::vector<double>& xs, std::vector<double>& ys) {
		// Testing function. For comparison only.
		// Fit centerPointX and centerPointY, which is exactly the floating-point result.
		// Without any local sampling result.

		// Prepare Parameters
		int rootsNum = xs.size();
		int paramsNum = highestDegree + 1;
	 	// Variables
		double xi, yi, ei, vari;
		double xij;
		gsl_matrix *X, *cov;
		gsl_vector *y, *w, *c;
		// Allocate
		// X, usually xi, xi^2, xi^3... (Expansion polynomial)
		X = gsl_matrix_alloc(rootsNum, paramsNum);
		// y, the observation of y value
		y = gsl_vector_alloc(rootsNum);
		// weight, usually the 1 / (err^2)
		w = gsl_vector_alloc(rootsNum);
		// store the results of polynomial parameters
		c = gsl_vector_alloc(paramsNum);
		cov = gsl_matrix_alloc(paramsNum, paramsNum);
		// Feed data into X, y, and w
		for (int i = 0; i < rootsNum; i++) {
			// Set X
			xi = xs[i];
			std::vector<double> terms_xi = {1.0, xi};
			terms_xi.resize(paramsNum);
			for (int j = 0; j < paramsNum; j++) {
				if (j > 1) {
					terms_xi[j] = terms_xi[j-1]*(2*xi)-terms_xi[j-2];
				}
				xij = terms_xi[j];
				gsl_matrix_set(X, i, j, xij);
			}
			// Set y
			yi = ys[i];
			gsl_vector_set(y, i, yi);
			gsl_vector_set(w, i, 1.0);
		}
		// Fitting using GSL's method
		gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(rootsNum, paramsNum);
		status = gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
		gsl_multifit_linear_free(work);
		if (status == 0) {
			// Store params from C to C++
			poly.coef.resize(paramsNum);
			for (int i = 0; i < paramsNum; i++) {
				poly.coef[i] = gsl_vector_get(c, i);
			}
			// Store params from C to C++
			xCov.resize(paramsNum);
			for (int i = 0; i < paramsNum; i++) {
				xCov[i].resize(paramsNum);
				for (int j = 0; j < paramsNum; j++) {
					xCov[i][j] = gsl_matrix_get(cov, i, j);
				}
			}
		}
		// Free space
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);
		return poly.coef;
	}
	int _fitVectors(const std::vector<double>& xList, const std::vector<double> yList) {
		// Testing function. For comparison only.
		// Prepare Parameters
		int rootsNum = xList.size();
		if (xList.size() != yList.size())
		{
			std::printf("Error in _fitVectors(), the length of xList and yList must be identical.\n");
			return -1;
		}
		int paramsNum = highestDegree + 1;
	 	// Variables
		double xi, yi, ei;
		double xij;
		gsl_matrix *X, *cov;
		gsl_vector *y, *w, *c;
		// Allocate
		// X, usually xi, xi^2, xi^3... (Expansion polynomial)
		X = gsl_matrix_alloc(rootsNum, paramsNum);
		// y, the observation of y value
		y = gsl_vector_alloc(rootsNum);
		// weight, usually the 1 / (err^2)
		w = gsl_vector_alloc(rootsNum);
		// store the results of polynomial parameters
		c = gsl_vector_alloc(paramsNum);
		cov = gsl_matrix_alloc(paramsNum, paramsNum);
		// Feed data into X, y, and w
		for (int i = 0; i < rootsNum; i++) {
			// Set X
			xi = xList[i];
			for (int j = 0; j < paramsNum; j++) {
				xij = poly.calcOnTerm(xi, j);
				gsl_matrix_set(X, i, j, xij);
			}
			// Set y
			yi = yList[i];
			gsl_vector_set(y, i, yi);
			// No weight information.
			gsl_vector_set(w, i, 1.0);
		}
		// Fitting using GSL's method
		gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(rootsNum, paramsNum);
		status = gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
		gsl_multifit_linear_free(work);
		if (status == 0) {
			// Store params from C to C++
			poly.coef.resize(paramsNum);
			for (int i = 0; i < paramsNum; i++) {
				poly.coef[i] = gsl_vector_get(c, i);
			}
			// Store params from C to C++
			xCov.resize(paramsNum);
			for (int i = 0; i < paramsNum; i++) {
				xCov[i].resize(paramsNum);
				for (int j = 0; j < paramsNum; j++) {
					xCov[i][j] = gsl_matrix_get(cov, i, j);
				}
			}
		}
		// Free space
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);
		return status;
	}
};

#endif