#ifndef LOCALSAMPLER_H
#define LOCALSAMPLER_H

#include <iostream>
#include <memory>
#include <cmath>
#include <random>
#include <iomanip>
#include <tuple>
#include <cassert>
#include <algorithm>
// #include <boost/math/statistics/bivariate_statistics.hpp>

#include "fpUtil.h"
#include "fpInterface.h"
#include "local-observation.h"

// #define LONGFP

#if defined (LONGFP)
typedef long double FPType;
#else
typedef double FPType;
#endif

class LocalSampler {
private:
    // Target function
    int funcIndex;
    std::unique_ptr<FloatingPointFunction> funcPtr;

    // Center Point
    double centerPointX;

    // Parameters on radius selection
    int insideSampleNum = 300;
    double minRadius = 1e-13;
    double maxRadius = 1e-1;
    int radiusNum = 100;

    // Filter extremely small radius, which contains only few x points.
    int radiusPointsThreshold = 100;

    // Parameters on intensive sampling with chosen radius
    int intensiveSampleNum = 50000;

    // Random device;
    std::random_device rd;
    std::mt19937_64 mt_generator_64;

public:
    LocalSampler(int index, double point) {
        funcIndex = index;
        funcPtr = std::make_unique<SimpleFunction>(funcIndex);
        centerPointX = point;
        // If use system random device as seed.
        mt_generator_64.seed(rd());
    }
public:
    // xMean, xRadius, yMean, yMin, yMax
    LocalObservation sampling() {
        // 1. Generate radius list
        std::vector<double> radiusList(radiusNum);
        double radiusLnOffset = (log(maxRadius)-log(minRadius)) / (radiusNum-1);
        for (int i = 0; i < radiusNum; i++) {
            double curRadius = exp(log(minRadius) + radiusLnOffset * i);
            radiusList[i] = curRadius;
            // std::cout << i << ' ' << std::scientific << curRadius << std::endl;
        }

        std::vector<double> varScoreList(radiusNum, 0.0);
        // 2. Sampling with each radius
        for (int i = 0; i < radiusNum; i++) {
            // 2.1 Define X list and Y list.
            std::vector<FPType> xList;
            std::vector<FPType> yList;

            // 2.2 Get sampling points for X list and Y list.
            double curRadius = radiusList[i];
            // std::cout << "Radius: " << curRadius << std::endl;

            int successCount = samplingWithRadius(curRadius, insideSampleNum, xList, yList);
            // -1 means the radius did not contain enough x points, filtered.
            if (successCount == -1) {
                continue;
            }

            // 2.3 Get basic statistic for X list and Y list.
            FPType xMin, xMax, yMin, yMax;
            xMin = xList.front();
            xMax = xList.back();
            // std::pair<FPType, FPType> xMM = getMinMax(xList);
            // xMin = xMM.first;
            // xMax = xMM.second;
            std::pair<FPType, FPType> yMM = getMinMax(yList);
            yMin = yMM.first;
            yMax = yMM.second;

            // 2.4.2 Diff-based measurement
            double curVarScore = calcAdjacentDiffScoreV3(xList, yList, xMin, xMax, yMin, yMax);
            // varScoreList[i] = curVarScore1 + curVarScore2 + curVarScore3;
            varScoreList[i] = curVarScore;

        }

        // 2.6 Select the best neighbour
        double bestRadius = 0;
        double bestCorr = 1.0;
        double bestVarScore = 0.0;

        // 2.6.3
        for (int i = 0; i < radiusNum; i++) {
            double curVarScore = varScoreList[i];
            if (curVarScore > bestVarScore) {
                bestVarScore = curVarScore;
                bestRadius = radiusList[i];
            }
        }

        // 3. Working on the chosen radius/neighbour.
        // 3.1 Show the chosen radius.
        // std::cout << std::scientific 
        //     << "\nCenter Point: " << centerPointX
        //     << "\nChosen Radius: " << bestRadius
        //     << std::endl;

        // 3.2 Intensive sampling in the target neighbour
        std::vector<FPType> intensiveXList;
        std::vector<FPType> intensiveYList;
        int intensiveSuccessCount = samplingWithRadius(bestRadius, intensiveSampleNum, intensiveXList, intensiveYList);
        (void)intensiveSuccessCount;
        // 3.3 Print result for analyzing.
        // auto [xMean, yMean, xyCov] = boost::math::statistics::means_and_covariance(intensiveXList, intensiveYList);
        double xMean = calcMean(intensiveXList);
        double yMean = calcMean(intensiveYList);
        double yVar = calcResidualVar(intensiveXList, intensiveYList);
        // std::pair<FPType, FPType> yMeanVar = calcMeanVar(intensiveYList);
        // double yMean = yMeanVar.first;
        // double yVar = yMeanVar.second;
        // auto[yMin, yMax]= getMinMax(intensiveYList);
        FPType yMin, yMax;
        std::pair<FPType, FPType> yMM = getMinMax(intensiveYList);
        yMin = yMM.first;
        yMax = yMM.second;
        // std::cout << std::scientific << std::setprecision(16)
        //     << "X Mean: " << xMean
        //     << "\n  Y Mean:   " << yMean
        //     << "   Y-Range:   " << yMax - yMin
        //     // << "\n  Y-Range Mean: " << (yMax+yMin)/2
        //     << std::endl;

        // 3.4 Get Floating Point Result for Center Point X.
        double centerPointY;
        funcPtr->call(centerPointX);
        centerPointY = funcPtr->getResult();

        LocalObservation res;
        res.xMean = xMean;
        res.yMean = yMean;
        res.yVar = yVar;
        res.radius = bestRadius;
        res.yMin = yMin;
        res.yMax = yMax;
        res.centerPointX = centerPointX;
        res.centerPointY = centerPointY;
        res.sampleNum = intensiveSampleNum;

        return res;
    }

private:
    double calcAdjacentDiffScoreV3(std::vector<FPType>& xList,
                        std::vector<FPType>& yList,
                        FPType xMin,
                        FPType xMax,
                        FPType yMin,
                        FPType yMax) {
        int len = xList.size();
        int len2 = yList.size();
        // assert(len == len2);
        assert(len > 1);

        FPType yTotalDiff = yMax - yMin;
        if (yTotalDiff == 0) return 0;

        FPType diff_threshold = yTotalDiff / len * 10; // 10 times larger than is should be
        FPType totalScore = 0.0, c = 0.0, offset, temp;
        for (int i = 0; i < len-1; i++) {
            if (fabs(yList[i+1]-yList[i]) < diff_threshold)
                continue;
            offset = fabs(yList[i+1]-yList[i]) - c;
            temp = totalScore + offset;
            c = (temp - totalScore) - offset;
            totalScore = temp;
        }
        totalScore = totalScore / (yMax - yMin);
        return totalScore;
    }

    FPType calcMean(const std::vector<FPType>& valList) {
        // Use Kahan Summation for sum.
        FPType sum = 0;
        FPType c = 0;
        FPType y, t;
        int len = valList.size();
        for (int i = 0; i < len; i++) {
            y = valList[i] - c;
            t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        FPType mean = sum / len;
        return mean;
    }

    std::pair<FPType, FPType> calcMeanVar(const std::vector<FPType>& valList) {
        // Get mean first.
        FPType mean = calcMean(valList);

        // Use mean to calculate variance.
        // Also apply Kahan Summation here.
        FPType sum2 = 0;
        FPType c = 0;
        FPType y, t;
        FPType df, dfsq;
        int len = valList.size();
        // Sum of (v[i] - mean)^2
        for (int i = 0; i < len; i++) {
            df = valList[i] - mean;
            dfsq = df * df;

            y = dfsq - c;
            t = sum2 + y;
            c = (t - sum2) - y;
            sum2 = t;
        }
        FPType var;
        if (len < 2)
            var = 0;
        else
            var = sum2 / (len - 1);

        return std::pair<FPType, FPType>(mean, var);
    }

    FPType calcResidualVar(const std::vector<FPType>& xList, const std::vector<FPType>& yList) {
        int ttlen = xList.size();
        int len[2] = {0,0};
        FPType sum[4] = {0,0,0,0}; // sumx0, sumx1, sumy0, sumy1
        FPType c[4] = {0,0,0,0};   // cx0, cx1, cy0, cy1
        FPType offset, temp;
        int ptr = 0;
        for (int i = 0; i < ttlen; i++) {
            if (xList[i] < centerPointX)
                ptr = 0;
            else
                ptr = 1;
            // count size
            len[ptr] += 1;
            // maintain sum of x
            offset = xList[i] - c[ptr];
            temp = sum[ptr] + offset;
            c[ptr] = (temp - sum[ptr]) - offset;
            sum[ptr] = temp;
            // maintain sum of y
            offset = yList[i] - c[ptr+2];
            temp = sum[ptr+2] + offset;
            c[ptr+2] = (temp - sum[ptr+2]) - offset;
            sum[ptr+2] = temp;
        }

        // two points to determine a linear function
        FPType x0 = sum[0] / len[0];
        FPType x1 = sum[1] / len[1];
        FPType y0 = sum[2] / len[0];
        FPType y1 = sum[3] / len[1];
        FPType k = (y1 - y0) / (x1 - x0);

        // calculate residual variance besides the linear explained.
        FPType varSum = 0;
        FPType varc = 0;
        FPType diff, var;
        for (int i = 0; i < ttlen; i++) {
            diff = yList[i] - (k * (xList[i] - x0) + y0);
            var = diff * diff;
            offset = var - varc;
            temp = varSum + offset;
            varc = (temp - varSum) - offset;
            varSum = temp;
        }
        if (ttlen < 2)
            return 0;
        else
            return varSum / (ttlen - 1);
    }

    std::pair<FPType, FPType> getMinMax(const std::vector<FPType>& valList) {
        if (valList.empty()) {
            return std::pair<FPType, FPType>(0,0);
        }
        FPType min = valList[0];
        FPType max = valList[0];
        int len = valList.size();
        for (int i = 1; i < len; i++) {
            if (valList[i] < min) {
                min = valList[i];
            }
            if (valList[i] > max) {
                max = valList[i];
            }
        }
        return std::pair<FPType, FPType>(min, max);
    }

    int samplingWithRadius(double radius, int sampleNum, std::vector<FPType>& xList, std::vector<FPType>& yList) {
        // Reserve size
        xList.reserve(sampleNum);
        yList.reserve(sampleNum);

        double fpLowerBound = centerPointX - fabs(centerPointX * radius);
        double fpUpperBound = centerPointX + fabs(centerPointX * radius);
        uint64_t i64LowerBound = fpUtil::doubleToI64(fpLowerBound);
        uint64_t i64UpperBound = fpUtil::doubleToI64(fpUpperBound);

        // swap for negative floating-point
        if (i64LowerBound > i64UpperBound) {
            uint64_t tmp = i64LowerBound;
            i64LowerBound = i64UpperBound;
            i64UpperBound = tmp;
        }

        // Filter extremely small radius, which contains only few x points.
        if (i64UpperBound - i64LowerBound < radiusPointsThreshold) {
            return -1;
        }
        
        // std::cout << "Lower bound: " << fpLowerBound << ' ' << i64LowerBound << std::endl;
        // std::cout << "Upper bound: " << fpUpperBound << ' ' << i64UpperBound << std::endl;
        // Now define in class member.
        // static std::random_device rd;
        // static std::mt19937_64 mt_generator_64(rd());
        // static std::mt19937_64 mt_generator_64(0xdeadbeef);
        std::uniform_int_distribution<uint64_t> uniDist(i64LowerBound, i64UpperBound);

        int successCount = 0;

        // Generate xList
        for (int i = 0; i < sampleNum; i++) {
            uint64_t i64x = uniDist(mt_generator_64);
            double fpx = fpUtil::i64ToDouble(i64x);
            xList.push_back(fpx);
        }
        // Sort xList
        sort(xList.begin(), xList.end());
        // Generate yList
        for (int i = 0; i < sampleNum; i++) {
            double fpx = xList[i];
            double fpy;
            funcPtr->call(fpx);
            if (!funcPtr->isSuccess()) {
                continue;
            }
            fpy = funcPtr->getResult();
            yList.push_back(fpy);
            successCount++;
        }
        return successCount;
    }
};

#endif
