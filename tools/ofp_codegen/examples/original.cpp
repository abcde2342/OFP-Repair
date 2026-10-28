#include <cmath>
#include <iomanip>
#include <iostream>

double cos_x2(double x) {
    return (1.0 - std::cos(x)) / (x * x);
}

int main() {
    std::cout << std::hexfloat;
    for (double x : {0.0, 1e-8, 1e-4, 0.01, 0.02}) {
        std::cout << x << '\t' << cos_x2(x) << '\n';
    }
}
