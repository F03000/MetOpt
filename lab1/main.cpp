#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

/// golden ratio for golden section method
const double GOLDEN_RATIO = (1 + sqrt(5)) / 2;

/// Saving number of iterations
size_t number_of_iterations;

/**
 * Function for research
 * @param x function argument
 * @return function value
 */
double func(double x) {
    return -3.0 * x * sin(0.75 * x) + exp(-2.0 * x);
}

/**
 * Dichotomy method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double dichotomy(double a, double b, double eps) {
    number_of_iterations = 0;
    double d = eps / 2;

    while (fabs(b - a) / 2 > eps) {
        number_of_iterations++;
        double x1 = (a + b) / 2 - d;
        double x2 = (a + b) / 2 + d;
        if (func(x1) <= func(x2)) {
            b = x2;
        } else {
            a = x1;
        }
    }

    return (a + b) / 2;
}

/**
 * Golden section method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double golden_section(double a, double b, double eps) {
    number_of_iterations = 0;

    while (fabs(b - a) / 2 > eps) {
        number_of_iterations++;
        double x1 = b - (b - a) / GOLDEN_RATIO;
        double x2 = a + (b - a) / GOLDEN_RATIO;
        if (func(x1) <= func(x2)) {
            b = x2;
        } else {
            a = x1;
        }
    }

    return (a + b) / 2;
}

double fib(int n) {
    return (pow((1 + sqrt(5)) / 2, n) - pow((1 - sqrt(5)) / 2, n)) / sqrt(5);
}

/**
 * Fibonacci method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param a0 left border
 * @param b0 right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double fibonacci(double a0, double b0, double eps) {
    std::vector<double> f(2);
    int n = 0;
    f[n++] = 1;
    f[n++] = 1;
    for (; f.back() <= (b0 - a0) / eps; ++n) {
        f.push_back(f.back() + f[n - 2]);
    }
    number_of_iterations = n - 2;

    double a = a0, b = b0;
    for (int k = 0; k < number_of_iterations; k++) {
        double x1 = a + (f[number_of_iterations - k - 1] * (b0 - a0)) / f.back();
        double x2 = a + (f[number_of_iterations - k] * (b0 - a0)) / f.back();
        if (func(x1) <= func(x2)) {
            b = x2;
        } else {
            a = x1;
        }
    }

    return (a + b) / 2;
}

/**
 * Parabolic method of finding min value
 * @param x1
 * @param x2
 * @param x3
 * @param prev_x
 * @param eps absolute accuracy
 * @return Min value in given range
 */
double parabolic(double x1, double x2, double x3, double prev_x, double eps) {
    std::cout << "a = " << x1 << ", b = " << x3 << std::endl;
    number_of_iterations++;
    double x1_res = func(x1), x2_res = func(x2), x3_res = func(x3);
    double a0 = x1_res, a1 = (x2_res - x1_res) / (x2 - x1), a2 =
            ((x3_res - x1_res) / (x3 - x1) - (x2_res - x1_res) / (x2 - x1)) / (x3 - x2);
    double x = (x1 + x2 - (a1 / a2)) / 2;
    double x_res = func(x);
    if (fabs(x - prev_x) <= eps) {
        return x;
    }
    if (x < x2) {
        if (x_res >= x2_res) {
            return parabolic(x, x2, x3, x, eps);
        } else {
            return parabolic(x1, x, x2, x, eps);
        }
    } else {
        if (x_res >= x2_res) {
            return parabolic(x1, x2, x, x, eps);
        } else {
            return parabolic(x2, x, x3, x, eps);
        }
    }
}

/**
 * Brent's method of finding min value
 * @require a < c && function unimodal in [a, c]
 * @param a left border
 * @param c right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double brent(double a, double c, double eps) {
    double x, w, v, x_res, w_res, v_res, d, e;
    x = v = w = (a + c) / 2;
    x_res = w_res = v_res = func(x);
    d = e = c - a;
    while (true) {
        double g = e, u = -1;
        e = d;
        if (x != v && x != w && w != v && x_res != v_res && x_res != w_res && v_res != w_res) {
            if (w > v) {
                double a0 = v_res, a1 = (x_res - v_res) / (x - v), a2 =
                        ((w_res - v_res) / (w - v) - (x_res - v_res) / (x - v)) / (w - x);
                u = (v + x - (a1 / a2)) / 2;
            } else {
                double a0 = w_res, a1 = (x_res - w_res) / (x - w), a2 =
                        ((v_res - w_res) / (v - w) - (x_res - w_res) / (x - w)) / (v - x);
                u = (w + x - (a1 / a2)) / 2;
            }
        }
        if (u - a > eps && u < c - eps && fabs(u - x) < g / 2) {
            d = fabs(u - x);
        } else {
            if (x < (c - a) / 2) {
                u = x + (c - x) / GOLDEN_RATIO;
                d = c - x;
            } else {
                u = x - (x - a) / GOLDEN_RATIO;
                d = x - a;
            }
        }
        if (std::abs(u - x) < eps) {
            return u;
//            u = x + ((u - x > 0) - (u - x < 0)) * eps;
        }
        double u_res = func(u);
        if (u_res <= x_res) {
            if (u >= x) {
                a = x;
            } else {
                c = x;
            }
            v = w;
            w = x;
            x = u;
            v_res = w_res;
            w_res = x_res;
            x_res = u_res;
        } else {
            if (u >= x) {
                c = u;
            } else {
                a = u;
            }
            if (u_res <= w_res || w == x) {
                v = w;
                w = u;
                v_res = w_res;
                w_res = u_res;
            } else if (u_res <= v_res || v == x || v == w) {
                v = u;
                v_res = u_res;
            }
        }
    }
}

/**
 * Console logger function
 * @param name: name of method
 * @param eps: accuracy
 * @param res: result
 */
void log(const std::string& name, double eps, double res) {
    std::cout << std::setprecision(8);
    std::cout << std::fixed;
    std::cout << "Method used:          " << name << std::endl;
    std::cout << "Absolute error        " << eps << std::endl;
    std::cout << "Result:               " << res << std::endl;
    std::cout << "Number of iterations: " << number_of_iterations << std::endl;
    std::cout << std::endl;
}

int main() {
    double eps = 10e-8;
    log("Dichotomy", eps, dichotomy(0, M_2_PI, eps));
    log("Golden section", eps, golden_section(0, M_2_PI, eps));
    log("Fibonacci", eps, fibonacci(0, M_2_PI, eps));
//    log("Parabolic", eps, parabolic());
//    log("Combined Brent", eps, brent(0, M_2_PI, eps));
    return 0;
}