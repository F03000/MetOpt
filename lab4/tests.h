#ifndef METOPT_TESTS_H
#define METOPT_TESTS_H

#include <iostream>
#include "methods.h"

void test01() {
    /*auto func2 = [](const vector_ &x) {
        return x[0] * x[0]  + x[1] * x[1] - 1.2 * x[0] * x[1];
    };
    auto grad2 = [](const vector_ &x) {
        return vector_ {2 * x[0] - 1.2 * x[1], 2 * x[1] - 1.2 * x[0]};
    };
    auto gess2 = [](const vector_ &x) {
        return matrix_ {vector_ {2, -1.2},
                        vector_ {-1.2, 2}};
    };
    extended_function f2 = extended_function(func2, grad2, gess2);
    vector_ x0_2 = vector_{4, 1};
    double eps = 10e-2;*/

    auto func2 = [](const vector_ &x) {
        return 100 * pow(x[0], 4) - 200 * pow(x[0], 2) * x[1]
               + pow(x[0], 2) - 2 * x[0] + 100 * pow(x[1], 2) + 1;
    };
    auto grad2 = [](const vector_ &x) {
        return vector_ {400 * pow(x[0], 3) - 400 * x[0] * x[1] + 2 * x[0] - 2, 200 * x[1] - 200 * pow(x[0], 2)};
    };
    auto gess2 = [](const vector_ &x) {
        return matrix_ {vector_ {1200 * pow(x[0], 2) - 400 * x[1] + 2, -400 * x[0]},
                        vector_ {-400 * x[0], 200.0}};
    };
    extended_function f2 = extended_function(func2, grad2, gess2);
    vector_ x0_2 = vector_{-1.2, 1};
    double eps = 10e-5;

    vector_ x1 = classic_newton(f2, x0_2, eps, 2);
    std::cout << std::endl << "classic: ";
    for (double i : x1) {
        std::cout << i << " ";
    }

    vector_ x2 = search_newton(f2, x0_2, eps, 2);
    std::cout << std::endl << "search: ";
    for (double i : x2) {
        std::cout << i << " ";
    }

    vector_ x3 = descent_newton(f2, x0_2, eps, 2);
    std::cout << std::endl << "descent: ";
    for (double i : x3) {
        std::cout << i << " ";
    }


//    vector_ x4 = dfp(f2, x0_2, eps);
//    vector_ x5 = pauel(f2, x0_2, eps);
}

/*// f(x) = x1*x2*x3
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
}*/

#endif //METOPT_TESTS_H
