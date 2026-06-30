#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <cstdio>
#include <vector>
#include <cassert>

#include "local-sampler.h"
#include "local-observation.h"
#include "interval-sampler.h"
#include "polynomial.h"
#include "serializer.h"

#define VERBOSE_ON_EACH_EVAL false

class Observation {
public:
	// Target function
	int funcIndex;

	// Interval under observe
	double leftBound;
	double rightBound;

	// Interval roots
	std::vector<double> roots;
	int rootsNum = 10000;

	// Observation
	std::vector<LocalObservation> rootsInfo;
public:
	Observation(int index, double lbound, double rbound) {
		funcIndex = index;
		leftBound = lbound;
		rightBound = rbound;
	}
	Observation(int index, double lbound, double rbound, int num) : Observation(index, lbound, rbound) {
		rootsNum = num;
	}
public:
	void init() {
		// 1. Calculate Chebyshev roots.
		IntervalSampler intervalSampler(funcIndex, leftBound, rightBound, rootsNum);
		roots = intervalSampler.getActualRoots();

		// 2. Sample on each root.
		rootsInfo.resize(rootsNum);

		std::printf("\nGenerating Observation... ");
		std::fflush(stdout);
		#pragma omp parallel for shared(rootsInfo) schedule(static)
		for (int i = 0; i < rootsNum; i++) {
			double root = roots[i];
			LocalSampler localSampler(funcIndex, root);
			rootsInfo[i] = localSampler.sampling();
			// show progress
			// if (i % 100 == 0)
			// 	std::printf("%d / %d\n", i, rootsNum);
		}
		std::printf("Finished.\n");
	}
public:
	int getSize() const { return rootsInfo.size(); }
	LocalObservation getObservation(int index) const { return rootsInfo[index]; }
	std::vector<LocalObservation> getObservations() const { return rootsInfo; }
	void writeObservationToFile() const {
		Serializer ser;
		ser.writeObsvToFile(rootsInfo);
	}
	void loadObservationFromFile() {
		Serializer ser;
		ser.readObsvFromFile(rootsInfo);
		rootsNum = rootsInfo.size();
	}

	double evaluatePoly(const Polynomial& poly) {
		if (VERBOSE_ON_EACH_EVAL) {
			std::printf("----------Start Evaluation of a Polynomial----------\n");
			std::printf("Score | X | Y on Poly | Y Mean | Y Variance | Radius\n");
		}

		// Calculate Score
		double totalScore = 0.0;
		for (int i = 0; i < rootsNum; i++) {
			// Get observation on each root
			const LocalObservation& curInfo = rootsInfo[i];
			double x = curInfo.xMean;

			// Get y on poly
			double yPoly = poly.evaluate(x);

			// Evaluate score for yPoly with observation
			double score = evaluatePolyPoint(yPoly, curInfo);

			totalScore += score;

			if (VERBOSE_ON_EACH_EVAL) {
				// Show polynomial's performance on this root
				printf("%.4f %.15e %.15e %.15e %.15e %.3e\n", score, curInfo.centerPointX, yPoly, curInfo.yMean, curInfo.yVar, curInfo.radius);
			}
		}

		if (VERBOSE_ON_EACH_EVAL) {
			// Show the Polynomial
			std::printf("\nEvaluated Polynomial: \n");
			poly.show();

			// And Overall Score
			std::printf("Score on above Polynomial: %.4f\n", totalScore);
			std::printf("----------------------------------------------------");
		}

		return totalScore;
	}
	double evaluateFunctionPtr(const std::unique_ptr<FloatingPointFunction>& ptr) {
		// similar with evaluatePoly
		if (VERBOSE_ON_EACH_EVAL) {
			std::printf("----------Start Evaluation of a Function Pointer----------\n");
			std::printf("Score | X | Y on Poly | Y Mean | Y Variance | Radius\n");
		}
		double totalScore = 0.0;
		for (int i = 0; i < rootsNum; i++) {
			const LocalObservation& curInfo = rootsInfo[i];
			double x = curInfo.xMean;
			double yFunc = ptr->callAndGetResult(x);
			double score = evaluatePolyPoint(yFunc, curInfo);
			totalScore += score;
			if (VERBOSE_ON_EACH_EVAL) {
				printf("%.4f %.15e %.15e %.15e %.15e %.3e\n", score, curInfo.centerPointX, yFunc, curInfo.yMean, curInfo.yVar, curInfo.radius);
			}
		}
		if (VERBOSE_ON_EACH_EVAL) {
			std::printf("Score on This Function Pointer: %.4f\n", totalScore);
			std::printf("----------------------------------------------------");
		}
		return totalScore;
	}

private:
	double evaluatePolyPoint(double yPoly, const LocalObservation& rootObs) const {
		double yMean = rootObs.yMean;
		double yMin = rootObs.yMin;
		double yMax = rootObs.yMax;
		// We need gradient to guide the search,
		// So negative score is still useful. 
		// 
		// if (yPoly < yMin || yPoly > yMax) {
		//     return 0.0;
		// }
		if (yPoly > yMean) {
			return 1.0 - (yPoly - yMean) / (yMax - yMean);
		}
		if (yPoly < yMean) {
			return 1.0 - (yMean - yPoly) / (yMean - yMin);
		}
		return 1.0;
	}
};

#endif