#ifndef METOPT_METHODS_H
#define METOPT_METHODS_H

#include <vector>
#include <functional>
#include <mpif-sizeof.h>
#include "linear_algebra.h"
#include "golden_section.h"


// TODO: научиться вычислять градиент и гессиан в точке
class extended_function {
public:

    extended_function() {

    }

    const vector_ grad(const vector_ &x) {

    }

    const matrix_ gess(const vector_ &x) {

        for (int i = 0; i < x.size(); ++i) {
            for (int j = 0; j < x.size(); ++j) {

            }
        }
    }
};

vector_ newton(const vector_ &x0, const double &eps,
               const std::function<vector_(vector_)> &getP,
               const std::function<double(vector_, vector_)> &getA) {
    vector_ x = vector_(x0.size());
    std::copy(x0.begin(), x0.end(), x.begin());
    while (true) {
        vector_ p = getP(x);
        double a = getA(x, p);
        vector_ pa = p * a;
        if (module(pa) < eps) {
            return x;
        }
        x = x + pa;
    }
}

vector_ classic_newton(extended_function &f, const vector_ &x0, const double &eps) {
    auto getP = [&](const vector_ &x){
            return reversed(f.gess(x)) * f.grad(x) * -1;
    };
    auto getA = [&](const vector_ &x, const vector_ &p) {
        return 1;
    };
    return newton(x0, eps, getP, getA);
}

vector_ search_newton(extended_function &f, const vector_ &x0, const double &eps) {
    auto getP = [&](const vector_ &x){
        return reversed(f.gess(x)) * f.grad(x) * -1;
    };
    auto getA = [&](const vector_ &x, const vector_ &p) {
        auto sf = [&](const double a) {
            return x - p * a;
        };
        return golden_section(sf, -1000, 1000, eps);
    };
    return newton(x0, eps, getP, getA);
}

vector_ descent_newton(extended_function &f, const vector_ &x0, const double &eps) {
    auto getP = [&](const vector_ &x){
        vector_ q = reversed(f.gess(x)) * f.grad(x) * -1;
        vector_ grad = f.grad(x);
        return scalar(q, grad) < 0 ? q : grad * -1;
    };
    auto getA = [&](const vector_ &x, const vector_ &p) {
        auto sf = [&](const double a) {
            return x - p * a;
        };
        return golden_section(sf, -1000, 1000, eps);
    };
    return newton(x0, eps, getP, getA);
}

vector_ quazinewton(extended_function &f, const vector_ &x0, const double &eps,
                    const std::function<matrix_ (vector_)> &getG) {
    matrix_ prev_g = matrix_(x0.size(), vector_(x0.size(), 0));
    for (int i = 0; i < x0.size(); ++i) {
        prev_g[i][i] = 1;
    }
    vector_ p = f.grad(x0) * -1;
    vector_ x = vector_(x0.size()), prev_x = vector_(x0.size());
    auto sf = [&](const double a) {
        return x - p * a;
    };
    double a = golden_section(sf, -1000, 1000, eps);
    std::copy(x.begin(), x.end(), prev_x.begin());
    x = x0 + p * a;

    while (true) {
        if (module(x - prev_x) < eps) {
            return x;
        }
        matrix_ g = getG(x);
        p = g * f.grad(x) * -1;
        a = golden_section(sf, -1000, 1000, eps);
        vector_ pa = p * a;
        std::copy(x.begin(), x.end(), prev_x.begin());
        prev_g = g;
        x = x + pa;
    }
}

vector_ dfp() {
    auto getG = [&](matrix_ prev_g, vector_ dx, vector_ dw) {
        int n = (int)prev_g.size();
        matrix_ first = matrix_(n, vector_(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                first[i][j] = dx[i] * dx[j];
            }
        }
        matrix_ second = matrix_(n, vector_(n));
        vector_ gw = prev_g * dw;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                second[i][j] = gw[i] * dw[j];
            }
        }
        return prev_g + first * (-1 / scalar(dw, dx)) +
                second * transparent(prev_g) * (-1 / scalar(prev_g * dw, dw));
    };
    return quazinewton();
}

vector_ pauel() {

}

#endif //METOPT_METHODS_H
