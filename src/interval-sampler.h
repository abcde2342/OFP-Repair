#ifndef INTERVALSAMPLER_H
#define INTERVALSAMPLER_H

#include <iostream>
#include <memory>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <cassert>
#include <algorithm>

#include "fpUtil.h"
#include "fpInterface.h"

const long double myPi = 3.14159265358979323846264338327950288419716939937510L;

class IntervalSampler {
private:
    // Target function
    int funcIndex;
    std::unique_ptr<FloatingPointFunction> funcPtr;

    // Interval under repair
    double leftBound;
    double rightBound;

    // Parameters on interval sampling
    int rootsNum = 10000;

    // Roots store
    std::vector<long double> nativeRoots;
    std::vector<double> actualRoots;

public:
    IntervalSampler(int index, double lbound, double rbound) {
        funcIndex = index;
        funcPtr = std::make_unique<SimpleFunction>(funcIndex);
        leftBound = lbound;
        rightBound = rbound;

        // Calculate nativeRoots
        calcNativeChebyshevRoots();
        calcActualChebyshevRoots();
    }

    IntervalSampler(int index, double lbound, double rbound, int num) {
        funcIndex = index;
        funcPtr = std::make_unique<SimpleFunction>(funcIndex);
        leftBound = lbound;
        rightBound = rbound;
        rootsNum = num;

        // Calculate nativeRoots
        calcNativeChebyshevRoots();
        calcActualChebyshevRoots();
    }

public:
    std::vector<long double> getNativeRoots() {
        return nativeRoots;
    }

    std::vector<double> getActualRoots() {
        return actualRoots;
    }
private:
    void calcNativeChebyshevRoots() {
        // x_{k}=\cos \left({\frac {2k-1}{2n}}\pi \right),\quad k=1,\ldots ,n.
        assert(rootsNum > 0);

        // clear vector to avoid duplicate.
        clearNativeRoots();

        // Using the direct formula for [-1, 1]
        // Maybe we can try some better one?
        for (int k = 1; k <= rootsNum; k++) {
            long double fpk = k;
            long double fpn = rootsNum;
            long double rt = cosl((2*fpk-1)/(2*fpn)*myPi);
            nativeRoots.push_back(rt);
        }
        std::sort(nativeRoots.begin(), nativeRoots.end());
    }

    void calcActualChebyshevRoots() {
        // Only called after native roots
        assert(nativeRoots.size() > 0);

        long double c1 = ((long double)leftBound + (long double)rightBound) / 2;
        long double c2 = ((long double)rightBound - (long double)leftBound) / 2;

        // clear vector to avoid duplicate.
        clearActualRoots();

        for (int i = 0; i < nativeRoots.size(); i++) {
            long double curNativeRoot = nativeRoots[i];
            long double curActualRoot = c1 + c2 * curNativeRoot;
            actualRoots.push_back(curActualRoot);
        }
    }

    void clearNativeRoots() {
        std::vector<long double> emptyList;
        std::swap(nativeRoots, emptyList);
    }
    void clearActualRoots() {
        std::vector<double> emptyList;
        std::swap(actualRoots, emptyList);
    }
};

#endif