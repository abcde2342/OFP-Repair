#ifndef SERIALIZER_H
#define SERIALIZER_H

#include <cstdio>
#include <string>
#include <cstring>
#include <cinttypes>

#include "local-observation.h"
#include "fpUtil.h"
#include "polynomial.h"

class Serializer {
private:
    static const std::string observationPath;

    static const std::string polynomialPath;

    static const std::string funcIDPath;

    static const std::string validationPath;

    static const std::string evaluationPath;

    static const std::string infoLevelPath;

public:
    static void cleanInfoLevelFile() {
        std::FILE* f = std::fopen(infoLevelPath.c_str(), "w");
        std::fclose(f);
        return;
    }
    static void writeInfoLevelToFile(int infolevel) {
        std::FILE* f = std::fopen(infoLevelPath.c_str(), "aw");
        std::fprintf(f, "%d\n", infolevel);
        std::fclose(f);
        return;
    }

    static void cleanEvaluationFile() {
        std::FILE* f = std::fopen(evaluationPath.c_str(), "w");
        std::fclose(f);
        return;
    }
    static void writeEvaluationToFile(
        double decayed_relative_origin,
        double decayed_relative_patch,
        double decayed_absolute_origin,
        double decayed_absolute_patch,
        double stable_relative_origin,
        double stable_relative_patch,
        double stable_absolute_origin,
        double stable_absolute_patch
    ) {
        std::FILE* f = std::fopen(evaluationPath.c_str(), "aw");

        std::fprintf(f, "%.16e,", decayed_relative_origin);
        std::fprintf(f, "%.16e,", decayed_relative_patch);
        std::fprintf(f, "%.16e,", decayed_absolute_origin);
        std::fprintf(f, "%.16e,", decayed_absolute_patch);
        std::fprintf(f, "%.16e,", stable_relative_origin);
        std::fprintf(f, "%.16e,", stable_relative_patch);
        std::fprintf(f, "%.16e,", stable_absolute_origin);
        std::fprintf(f, "%.16e\n", stable_absolute_patch);

        std::fclose(f);
        return;
    }

    static void cleanValidationFile() {
        std::FILE* f = std::fopen(validationPath.c_str(), "w");
        std::fclose(f);
        return;
    }
    static void writeValidationToFile(bool valid) {
        std::FILE* f = std::fopen(validationPath.c_str(), "aw");
        std::fprintf(f, "%d\n", valid);
        std::fclose(f);
        return;
    }

    static void cleanFuncIDFile() {
        std::FILE* f = std::fopen(funcIDPath.c_str(), "w");
        std::fclose(f);
        return;
    }
    static void writeFuncIDToFile(std::string funcID) {
        std::FILE* f = std::fopen(funcIDPath.c_str(), "aw");
        std::fprintf(f, "%s\n", funcID.c_str());
        std::fclose(f);
        return;
    }

    static void writePolyToFile(
        const int mode,
        const double lbound,
        const double rbound, 
        const double c1,
        const double c2,
        const std::vector<double>& param)
    {

        int num = param.size();

        std::FILE* f = std::fopen(polynomialPath.c_str(), "w");

        uint64_t mode_i64 = (uint64_t)mode;
        uint64_t lb_i64 = fpUtil::doubleToI64(lbound);
        uint64_t rb_i64 = fpUtil::doubleToI64(rbound);
        uint64_t c1_i64 = fpUtil::doubleToI64(c1);
        uint64_t c2_i64 = fpUtil::doubleToI64(c2);

        std::fprintf(f, "%" PRIx64 "\n", mode_i64);
        std::fprintf(f, "%" PRIx64 "\n", lb_i64);
        std::fprintf(f, "%" PRIx64 "\n", rb_i64);
        std::fprintf(f, "%" PRIx64 "\n", c1_i64);
        std::fprintf(f, "%" PRIx64 "\n", c2_i64);

        for (int i = 0; i < num; i++) {
            uint64_t p_i64 = fpUtil::doubleToI64(param[i]);
            std::fprintf(f, "%" PRIx64 "\n", p_i64);
        }

        std::fclose(f);

        return;
    }

    static void readPolyFromFile(
        int& mode,
        double& lbound,
        double& rbound,
        double& c1,
        double& c2,
        std::vector<double>& param)
    {
        param.clear();

        std::FILE* f = std::fopen(polynomialPath.c_str(), "r");

        if (f == NULL) {
            std::printf("File does not exist: %s\n", polynomialPath.c_str());
            return;
        }

        const int bufferSize = 256;
        char buffer[bufferSize];

        int linecount = 0;
        while ( (std::fgets(buffer, bufferSize, f) != NULL) ) {
            uint64_t p_i64;
            std::sscanf(buffer, "%" PRIx64, &p_i64);

            if (linecount == 0) {
                mode = (int)p_i64;
            }
            else if (linecount == 1) {
                lbound = fpUtil::i64ToDouble(p_i64);
            }
            else if (linecount == 2) {
                rbound = fpUtil::i64ToDouble(p_i64);
            }
            else if (linecount == 3) {
                c1 = fpUtil::i64ToDouble(p_i64);
            }
            else if (linecount == 4) {
                c2 = fpUtil::i64ToDouble(p_i64);
            }
            else {
                param.push_back(fpUtil::i64ToDouble(p_i64));
            }
            linecount++;
        }
        return;
    }

    static void writeObsvToFile(const std::vector<LocalObservation>& infoList) {
        int num = infoList.size();

        std::FILE* f = std::fopen(observationPath.c_str(), "w");

        for (int i = 0; i < num; i++) {
            LocalObservation info = infoList[i];
            uint64_t xMean_i64          = fpUtil::doubleToI64(info.xMean);
            uint64_t yMean_i64          = fpUtil::doubleToI64(info.yMean);
            uint64_t yVar_i64           = fpUtil::doubleToI64(info.yVar);
            uint64_t radius_i64         = fpUtil::doubleToI64(info.radius);
            uint64_t yMin_i64           = fpUtil::doubleToI64(info.yMin);
            uint64_t yMax_i64           = fpUtil::doubleToI64(info.yMax);
            uint64_t centerPointX_i64   = fpUtil::doubleToI64(info.centerPointX);
            uint64_t centerPointY_i64   = fpUtil::doubleToI64(info.centerPointY);
            uint64_t sampleNum_i64      = fpUtil::doubleToI64(info.sampleNum);
            std::fprintf(f, "%" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 "\n",
                xMean_i64,
                yMean_i64,
                yVar_i64,
                radius_i64,
                yMin_i64,
                yMax_i64,
                centerPointX_i64,
                centerPointY_i64,
                sampleNum_i64
            );
        }

        std::fclose(f);

        return;
    }

    static void readObsvFromFile(std::vector<LocalObservation>& infoList) {
        infoList.clear();

        std::FILE* f = std::fopen(observationPath.c_str(), "r");

        if (f == NULL) {
            std::printf("File does not exist: %s\n", observationPath.c_str());
            return;
        }

        const int bufferSize = 256;
        char buffer[bufferSize];

        while ( (std::fgets(buffer, bufferSize, f) != NULL) ) {
            uint64_t xMean_i64;
            uint64_t yMean_i64;
            uint64_t yVar_i64;
            uint64_t radius_i64;
            uint64_t yMin_i64;
            uint64_t yMax_i64;
            uint64_t centerPointX_i64;
            uint64_t centerPointY_i64;
            uint64_t sampleNum_i64;

            std::sscanf(buffer, "%" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64,
                &xMean_i64,
                &yMean_i64,
                &yVar_i64,
                &radius_i64,
                &yMin_i64,
                &yMax_i64,
                &centerPointX_i64,
                &centerPointY_i64,
                &sampleNum_i64
            );

            LocalObservation info;
            info.xMean          = fpUtil::i64ToDouble(xMean_i64);
            info.yMean          = fpUtil::i64ToDouble(yMean_i64);
            info.yVar           = fpUtil::i64ToDouble(yVar_i64);
            info.radius         = fpUtil::i64ToDouble(radius_i64);
            info.yMin           = fpUtil::i64ToDouble(yMin_i64);
            info.yMax           = fpUtil::i64ToDouble(yMax_i64);
            info.centerPointX   = fpUtil::i64ToDouble(centerPointX_i64);
            info.centerPointY   = fpUtil::i64ToDouble(centerPointY_i64);
            info.sampleNum      = fpUtil::i64ToDouble(sampleNum_i64);

            infoList.push_back(info);
        }

        return;
    }
};

const std::string Serializer::observationPath = "data/observation.dat";
const std::string Serializer::polynomialPath  = "data/polynomial.dat";
const std::string Serializer::funcIDPath      = "data/funcID.dat";
const std::string Serializer::validationPath  = "data/validation.dat";
const std::string Serializer::evaluationPath  = "data/evaluation.dat";
const std::string Serializer::infoLevelPath   = "data/infolevel.dat";

#endif