#include "gsl-fit.h"

int main(int argc, char **argv) {

	// Do the sampling to get the observations
    int funcIndex = (argc >= 2) ? atoi(argv[1]) : 2;
    double center = (argc >= 3) ? atof(argv[2]) : 0;
    double lbound = center + ((argc >= 4) ? atof(argv[3]) : -1e-2);
    double rbound = center + ((argc >= 5) ? atof(argv[4]) : 1e-2);
	int rootsNum = (argc >= 6) ? atoi(argv[5]) : 1000;

    // Try different settings
    int highestDegree = (argc >= 7) ? atoi(argv[6]) : 6;
    int weighted = (argc >= 8) ? atoi(argv[7]) : 1;
    int polymode = MODE_CHEBY; // MODE_POWER(0) or MODE_CHEBY(1)
    
    const char *fitmode = getenv("FITMODE");
    // const char *fitmode = "CHEBY_ONLY";
    // Possible Choice:
    // CHEBY_WEIGHTED_OBSV
    // CHEBY_OBSV
    // CHEBY_ONLY
    // NONE
    if (fitmode)
        printf("[OK] Read env. var. FITMODE=\"%s\"\n", fitmode);
    
    const char *polymode_env = getenv("POLYMODE");
    if (polymode_env) {
        if (!strcmp(polymode_env, "MODE_POWER")) polymode = MODE_POWER;
        if (!strcmp(polymode_env, "MODE_CHEBY")) polymode = MODE_CHEBY;
        printf("[OK] Read env. var. POLYMODE=\"%s\"\n", polymode_env);
    }

    Polynomial poly;

	double score;
	Observation observation(funcIndex, lbound, rbound, rootsNum);
    observation.init();

    // Write function index to file
    Serializer ser;
    // ser.writeFuncIDToFile(funcIndex);

    // Write observation to file
    observation.writeObservationToFile();

    // Use GSL fit to get initial fit results
    GSLFitModule gslfit;

    if (fitmode == NULL || !strcmp(fitmode, "CHEBY_WEIGHTED_OBSV")) {
        // Test Weighted
        // 1. Using GSL's fitting method to generate polynomial
        printf("GSL's Weighted Fitting...\n");
        if (fitmode) {
            weighted = 1;
        }
        gslfit.init(highestDegree, weighted, polymode, lbound, rbound);
        gslfit.fitObservations(observation);
        poly = gslfit.getPoly();
        poly.checkacc0();
        printf("Generated Polynomial: \n");
        poly.show();
        poly.writePolyToFile();
        // 2. Evaluation the poly on observations
        score = observation.evaluatePoly(poly);
        printf("\nScore on above Poly: %.4f\n", score/rootsNum);
        if (score > 0) {
            poly.writeValidToFile(1);
            printf("Valid.\n");
        }
        else {
            poly.writeValidToFile(0);
            printf("Invalid.\n");
        }
        printf("#####################################################\n");
    }

    if (fitmode && !strcmp(fitmode, "CHEBY_OBSV")) {
        // Test Un-Weighted
        // printf("GSL's Un-Weighted Fitting...\n");
        gslfit.init(highestDegree, 0, polymode, lbound, rbound);
        gslfit.fitObservations(observation);
        poly = gslfit.getPoly();
        // printf("Generated Polynomial: \n");
        poly.show();
        poly.writePolyToFile();
        score = observation.evaluatePoly(poly);
        // printf("\nScore on above Poly: %.4f\n", score/rootsNum);
        // printf("#####################################################\n");
    }
    
    // if (fitmode && !strcmp(fitmode, "CHEBY_ONLY")) {
    //     // Test with Chebyshev roots, without local sampling 
    //     // printf("Chebyshev roots, without local sampling...\n");
    //     gslfit.init(highestDegree, 0, polymode, lbound, rbound);
    //     gslfit._fitCenterPointsOnly(observation);
    //     poly = gslfit.getPoly();
    //     // printf("Generated Polynomial: \n");
    //     poly.show();
    //     poly.writePolyToFile();
    //     score = observation.evaluatePoly(poly);
    //     // printf("\nScore on above Poly: %.4f\n", score/rootsNum);
    //     // printf("#####################################################\n");
    // }
 
    if (fitmode && !strcmp(fitmode, "NONE")) {
        // Test without Chebyshev roots, without local sampling
        // manually generate data here
        // 1. random utilities
        // printf("Random roots, without local sampling...\n");
        std::random_device rd;
        std::mt19937_64 mt_generator_64;
        mt_generator_64.seed(rd());
        std::uniform_real_distribution<double> uniDist(lbound, rbound);
        // 2. function pointer
        std::unique_ptr<FloatingPointFunction> funcPtr = std::make_unique<SimpleFunction>(funcIndex);
        // 3. generate data
        std::vector<double> xList;
        std::vector<double> yList;
        for (int i = 0; i < rootsNum; i++) {
            double x, y;
            x = uniDist(mt_generator_64);
            funcPtr->call(x);
            if (!funcPtr->isSuccess()) {
                std::cout << "Error Input: " << x << std::endl;
                std::cin.get();
                continue;
            }
            y = funcPtr->getResult();
            xList.push_back(x);
            yList.push_back(y);
        }
        // 4. gsl fit
        gslfit.init(highestDegree, 0, polymode, lbound, rbound);
        gslfit._fitVectors(xList, yList);
        poly = gslfit.getPoly();
        // printf("Generated Polynomial: \n");
        poly.show();
        poly.writePolyToFile();
        score = observation.evaluatePoly(poly);
        // printf("\nScore on above Poly: %.4f\n", score/rootsNum);
        // printf("#####################################################\n");
    }

	return 0;
}
