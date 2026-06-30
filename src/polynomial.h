#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cmath>
#include <cstdio>
#include "serializer.h"
#include <vector>
#include "gsl-fit2.h"

#define MODE_POWER 0
#define MODE_CHEBY 1
#define MODE_CHEBYACC 2

class Polynomial {
public:
	int mode; // 0 for power series, 1 for chebyshev series

	std::vector<double> coef; // coefficients

	/* For range reduction.
	 * x in [l, r] => u in [-1, 1]
	 * u = (2.0*x-l-r) / (r-l)
	 * out of bound check: if ((x-a)*(x-b)>0.0) throw out of bound error
	 */
	double lbound, rbound;

	double c1, c2;

	Polynomial() {
		mode = 0;
		lbound = -1;
		rbound = 1;
		c1 = 0.0;
		c2 = 0.0;
	}

	Polynomial(int md, double lb, double rb) {
		mode = md;
		lbound = lb;
		rbound = rb;
		c1 = 0.0;
		c2 = 0.0;
	}

	void writeValidToFile(bool valid) const {
		Serializer ser;
		ser.writeValidationToFile(valid);
	}

	void writePolyToFile() const {
		Serializer ser;
		ser.writePolyToFile(mode, lbound, rbound, c1, c2, coef);
	}
	void loadPolyFromFile() {
		Serializer ser;
		ser.readPolyFromFile(mode, lbound, rbound, c1, c2, coef);
	}

	double reduction(double x) const {
		// long double lx = x;
		// long double llb = lbound;
		// long double lrb = rbound;
		// long double lu = (2.0*x - llb - lrb) / (lrb - llb);

		double sum = 2*x, c = 0.0, temp, u;
		temp = sum - lbound;
		if (std::abs(sum) >= std::abs(lbound))
			c = c + ((sum-temp)-lbound);
		else
			c = c + ((-lbound-temp)+sum);
		sum = temp;
		temp = sum - rbound;
		if (std::abs(sum) >= std::abs(rbound))
			c = c + ((sum-temp)-rbound);
		else
			c = c + ((-rbound-temp)+sum);
		sum = temp;
		sum = sum + c;
		u = sum / (rbound - lbound);

		return u;
	}

	double evaluate(double x) const {
		if (mode == MODE_CHEBYACC)
			return chebAccEval(x);
		else if (mode == MODE_CHEBY)
			return chebEvaluate(x);
		return powerEvaluate(x);
	}

	double Clenshaw(double x) const {
		double d = 0.0;
		double dd = 0.0;

		int nums = coef.size();
		for (int i = nums-1; i>=1; i--) {
			double tmp = d;
			d = 2*x*d - dd + coef[i];
			dd = tmp;
		}
		d = x*d - dd + coef[0];
		return d;
	}

	double chebEvaluate(double x) const {
		int nums = coef.size();

		double d = 0.0;
		double dd = 0.0;

		double u = reduction(x);
		double u2 = 2.0*u;

		// Clenshaw method
		for (int i = nums-1; i>=1; i--) {
			double tmp = d;
			d = u2*d - dd + coef[i];
			dd = tmp;
		}
		d = u*d - dd + coef[0];
		return d;
	}

	double powerEvaluate(double x) const {
		int nums = coef.size();

		double u = reduction(x);
		double y;
		// Expansion series
		// https://en.wikipedia.org/wiki/Horner%27s_method
		// lower-order coef first:
		// coef[0] + coef[1]*x + coef[2]*x^2 + coef[3]^3
		y = coef[nums-1];
		for (int i = nums-2; i >= 0; i--) {
			y = y * u + coef[i];
		}
		return y;
	}

	std::vector<long double> vandermonde(double x, int cols) {
		if (mode == MODE_CHEBY || mode == MODE_CHEBYACC)
			return chebVander(x, cols);
		return powerVander(x, cols);
	}

	std::vector<long double> powerVander(double x, int cols) {
		long double u = reduction(x);
		std::vector<long double> vander(cols);
		if (cols == 0)
			return vander;
		vander[0] = 1.0;
		for (int i = 1; i < cols; i++) {
			vander[i] = vander[i-1] * u;
		}
		return vander;
	}

	std::vector<long double> chebVander(double x, int cols) {
		long double u = reduction(x);
		std::vector<long double> vander(cols);
		if (cols == 0)
			return vander;
		vander[0] = 1.0;
		if (cols == 1)
			return vander;
		vander[1] = u;
		for (int i = 2; i < cols; i++) {
			vander[i] = vander[i-1] * (2 * u) - vander[i-2];
		}
		return vander;
	}

private:
	double chebAccEval(double x) const {
		if (fabs(x) < 1.4901161193847656e-08) {
			double val = x*(c1 + x*c2);
			return val;
		}
		else
			return chebEvaluate(x);
	}
public:
	double calcOnTerm(double x, int termNum) const {
		double u = reduction(x);
		// Expansion series
		return pow(u, termNum);
	}

	void checkacc0() {
		printf("PolyCheck: mode: %d, l: %f, r: %f\n", mode, lbound, rbound);
		if (mode == MODE_CHEBY) {
			if (lbound < 0 && rbound > 0)
				// printf("PolyCheck: %.15e\n", fabs(evaluate(0.0)));
				// getchar();
				if (fabs(evaluate(0.0)) < 1e-13) {
					mode = MODE_CHEBYACC;
					double l = std::max(-1e-4, lbound);
					double r = std::min(1e-4, rbound);
					int nums = 1000;
					double gap = (r-l) / nums;
					std::vector<double> xs, ys;
					for (int i = 0; i < nums; i++) {
						double xi = l + i*gap;
						double yi = evaluate(xi);
						xs.push_back(xi);
						ys.push_back(yi);
					}
					std::pair<double,double> c12 = fitPower2(xs, ys);
					c1 = c12.first;
					c2 = c12.second;
				}
		}
	}

	void show() const {
		// Show the type of polynomial
		if (mode == MODE_POWER) 
			printf("Polynomial Type: Expansion: \n");
		if (mode == MODE_CHEBY)
			printf("Polynomial Type: Chebyshev: \n");
		if (mode == MODE_CHEBYACC) {
			printf("Polynomial Type: Chebyshev Acc: \n");
			printf("c1: %.5e\nc2: %.5e\n", c1, c2);
		}

		// Show the coefficients of polynomial
		std::printf("[\n");
		for (int i = 0; i < coef.size(); ++i)
		{
			printf("  %.16e,\n", coef[i]);
		}
		std::printf("]\n");
	}
};

#endif