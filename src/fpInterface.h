#ifndef FPINTERFACE_RP_H
#define FPINTERFACE_RP_H

#include "../inputs/target.h"

class FloatingPointFunction {
public:
    virtual void call(double x) = 0;
    virtual double getResult() = 0;
    virtual double callAndGetResult(double x) = 0;
    virtual bool isSuccess() = 0;
    virtual ~FloatingPointFunction() {}

protected:
    double in, out;
    int status;
};

class SimpleFunction : public FloatingPointFunction {
public:
    SimpleFunction(int index) {
        if (index < 0 || index >= simpleFuncList.size()) {
            std::cout << "Invalid index in [SimpleFunction]: " << index << '\n';
            funcRef = simpleFuncList[0];
            return;
        }
        funcRef = simpleFuncList[index];
    }
    void call(double x) {
        in = x;
        out = funcRef(x);
        if (std::isnan(out) || std::isinf(out)) {
            status = -1;
        }
        else {
            status = 0;
        }
    }
    double callAndGetResult(double x) {
        call(x);
        return out;
    }

    double getResult() { return out; }
    bool isSuccess() { return (status == 0); }

private:
    std::function<double(double)> funcRef;
};

#endif
