#include "tests.h"

int main() {
    extended_function f;
    vector_ x0 = vector_(2, 9);
    double eps = 10e-7;

    vector_ x1 = classic_newton(f, x0, eps);
    vector_ x2 = search_newton(f, x0, eps);
    vector_ x3 = descent_newton(f, x0, eps);

    vector_ x4 = dfp(f, x0, eps);
    vector_ x5 = pauel(f, x0, eps);
}