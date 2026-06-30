#ifndef LOCALOBSERVATION_H
#define LOCALOBSERVATION_H

struct LocalObservation {
    double xMean;
    double yMean;
    double yVar;
    double radius;
    double yMin;
    double yMax;
    // Only for debug
    // Notice, centerX not equals xMean. 
    double centerPointX;
    double centerPointY;
    double sampleNum;
    bool valid;

	bool operator==(const LocalObservation& b) const {
    	return (
        (xMean == b.xMean) &&
        (yMean == b.yMean) &&
        (yVar == b.yVar) &&
        (radius == b.radius) &&
        (yMin == b.yMin) &&
        (yMax == b.yMax) &&
        (centerPointX == b.centerPointX) &&
        (centerPointY == b.centerPointY) &&
        (sampleNum == b.sampleNum) &&
        (valid == b.valid)
    	);
	}
};

#endif