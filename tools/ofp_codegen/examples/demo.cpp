#include <cmath>
#include <iomanip>
#include <iostream>

#include "cos_x2.hpp"

double original(double x) {
    return (1.0 - std::cos(x)) / (x * x);
}

double repaired(double x) {
    return ofp_codegen::cos_x2_repaired(x, &original);
}

int main() {
    std::cout << std::setprecision(17);
    for (double x : {0.0, 1e-8, 1e-4, 0.02}) {
        std::cout << "x=" << x << " original=" << original(x)
                  << " repaired=" << repaired(x) << '\n';
    }
}
