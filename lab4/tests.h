#ifndef METOPT_TESTS_H
#define METOPT_TESTS_H

#include "methods.h"

// f(x) = x1^3 + x2^3 + x3^3 - 5
void test01() {
    auto func = [](const vector_ &x) {
        return x[0] * x[0] * x[0] + x[1] * x[1] * x[1] + x[2] * x[2] * x[2] - 5;
    };
    auto grad = [](const vector_ &x) {
        return vector_ {3 * x[0] * x[0], 3 * x[0] * x[0], 3 * x[0] * x[0]};
    };
    auto gess = [](const vector_ &x) {
        return matrix_ {vector_ {6 * x[0], 0, 0},
                        vector_ {0, 6 * x[1], 0},
                        vector_ {0, 0, 6 * x[2]}};
    };
    extended_function f = extended_function(func, grad, gess);
    vector_ x0 = vector_(3, 0);
    double eps = 10e-7;

    vector_ x1 = classic_newton(f, x0, eps);
    vector_ x2 = search_newton(f, x0, eps);
    vector_ x3 = descent_newton(f, x0, eps);

    vector_ x4 = dfp(f, x0, eps);
    vector_ x5 = pauel(f, x0, eps);
}

// f(x) = x1*x2*x3
void test02() {
    auto func = [](const vector_ &x) {
        return x[0] * x[1] * x[2];
    };
    auto grad = [](const vector_ &x) {
        return vector_ {x[2] * x[3], x[1] * x[3], x[1] * x[2]};
    };
    auto gess = [](const vector_ &x) {
        return matrix_ {vector_ {0, x[3], x[2]},
                        vector_ {x[3], 0, x[1]},
                        vector_ {x[2], x[1], 0}};
    };
    extended_function f = extended_function(func, grad, gess);
    vector_ x0 = vector_(3, 0);
    double eps = 10e-7;

    vector_ x1 = classic_newton(f, x0, eps);
    vector_ x2 = search_newton(f, x0, eps);
    vector_ x3 = descent_newton(f, x0, eps);

    vector_ x4 = dfp(f, x0, eps);
    vector_ x5 = pauel(f, x0, eps);
}

#endif //METOPT_TESTS_H
