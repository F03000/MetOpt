#ifndef METOPT_METHODS_H
#define METOPT_METHODS_H

#include <vector>
#include <functional>
#include <mpif-sizeof.h>
#include "linear_algebra.h"
#include "golden_section.h"

class extended_function {
private:
    std::function<matrix_(vector_)> _gess;
    std::function<vector_(vector_)> _grad;
    std::function<double(vector_)> _func;
public:
    extended_function(std::function<double(vector_)> &func,
                      std::function<vector_(vector_)> &grad,
                      std::function<matrix_(vector_)> &gess) {
        _func = func;
        _gess = gess;
        _grad = grad;
    }

    vector_ grad(const vector_ &x) {
        return _grad(x);
    }

    matrix_ gess(const vector_ &x) {
        return _gess(x);
    }

    double func(const vector_ &x) {
        return _func(x);
    }
};

// TODO: сделать вывод в файл в методах newton и quazinewton
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
            return x + p * a;
        };
        return golden_section(sf, eps);
    };
    return newton(x0, eps, getP, getA);
}

vector_ descent_newton(extended_function &f, const vector_ &x0, const double &eps) {
    // TODO: поиск направления спуска через СЛАУ
    auto getP = [&](const vector_ &x){
        vector_ q = reversed(f.gess(x)) * f.grad(x) * -1;
        vector_ grad = f.grad(x);
        return scalar(q, grad) < 0 ? q : grad * -1;
    };
    auto getA = [&](const vector_ &x, const vector_ &p) {
        auto sf = [&](const double a) {
            return x - p * a;
        };
        return golden_section(sf, eps);
    };
    return newton(x0, eps, getP, getA);
}

vector_ quazinewton(extended_function &f, const vector_ &x0, const double &eps,
                    const std::function<matrix_ (matrix_, vector_, vector_)> &getG) {
    matrix_ prev_g = matrix_(x0.size(), vector_(x0.size(), 0));
    for (int i = 0; i < x0.size(); ++i) {
        prev_g[i][i] = 1;
    }
    vector_ prev_w = f.grad(x0) * -1;
    vector_ p = prev_w, x = x0;
    auto sf = [&](const double a) {
        return x - p * a;
    };
    double a = golden_section(sf, eps);
    vector_ prev_x = x;
    x = x0 + p * a;

    while (true) {
        vector_ dx = x - prev_x;
        if (module(x - prev_x) < eps) {
            return x;
        }
        vector_ w = f.grad(x);
        vector_ dw = w - prev_w;
        matrix_ g = getG(prev_g, dx, dw);
        p = g * f.grad(x) * -1;
        a = golden_section(sf, eps);
        prev_g = g;
        prev_x = x;
        prev_w = w;
        x = x + p * a;
    }
}

vector_ dfp(extended_function &f, const vector_ &x0, const double eps) {
    auto getG = [&](const matrix_ &prev_g, const vector_ &dx, const vector_ &dw) {
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
    return quazinewton(f, x0, eps, getG);
}

vector_ pauel(extended_function &f, const vector_ &x0, const double eps) {
    auto getG = [&](const matrix_ &prev_g, const vector_ &dx, const vector_ &dw) {
        int n = (int)prev_g.size();
        vector_ y = dx + prev_g * dw;
        matrix_ m = matrix_(n, vector_(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                m[i][j] = y[i] * y[j];
            }
        }
        return prev_g + m * (-1 / scalar(dw, y));
    };
    return quazinewton(f, x0, eps, getG);
}

#endif //METOPT_METHODS_H
