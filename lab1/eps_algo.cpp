#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <functional>
#include <fstream>

typedef const std::function<double(double)>& func;

/// golden ratio for golden section method
const double GOLDEN_RATIO = (1 + sqrt(5)) / 2;

/// Saving number of iterations
size_t number_of_iterations;

/**
 * Dichotomy method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param f function for research
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double dichotomy(func f, double a, double b, double eps) {
    number_of_iterations = 0;
    double d = eps / 2;
    while (fabs(b - a) / 2 > eps) {
        number_of_iterations++;
        double x1 = (a + b) / 2 - d;
        double x2 = (a + b) / 2 + d;
        if (f(x1) <= f(x2)) {
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
 * @param f function for research
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double golden_section(func f, double a, double b, double eps) {
    number_of_iterations = 0;
    while (fabs(b - a) / 2 > eps) {
        number_of_iterations++;
        double x1 = b - (b - a) / GOLDEN_RATIO;
        double x2 = a + (b - a) / GOLDEN_RATIO;
        if (f(x1) <= f(x2)) {
            b = x2;
        } else {
            a = x1;
        }
    }
    return (a + b) / 2;
}

/**
 * Fibonacci method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param f function for research
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double fibonacci(func f, double a0, double b0, double eps) {
    std::vector<double> fib(2);
    int n = 0;
    fib[n++] = 1;
    fib[n++] = 1;
    for (; fib.back() <= (b0 - a0) / eps; ++n) {
        fib.push_back(fib.back() + fib[n - 2]);
    }
    number_of_iterations = n - 2;

    double a = a0, b = b0;
    for (int k = 0; k < number_of_iterations; k++) {
        double x1 = a + (fib[number_of_iterations - k - 1] * (b0 - a0)) / fib.back();
        double x2 = a + (fib[number_of_iterations - k] * (b0 - a0)) / fib.back();
        if (f(x1) <= f(x2)) {
            b = x2;
        } else {
            a = x1;
        }
    }
    return (a + b) / 2;
}

/**
 * Parabolic method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param f function for research
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double parabolic(func f, double a, double b, double eps) {
    double prev_x = a, x1 = a, x2 = (a + b) / 2, x3 = b;
    double f_x1 = f(x1), f_x2 = f(x2), f_x3 = f(x3);
    number_of_iterations = 0;
    while (true) {
        number_of_iterations++;
        double a0 = f_x1, a1 = (f_x2 - f_x1) / (x2 - x1), a2 =
                ((f_x3 - f_x1) / (x3 - x1) - (f_x2 - f_x1) / (x2 - x1)) / (x3 - x2);
        double x = (x1 + x2 - (a1 / a2)) / 2;
        double f_x = f(x);
        if (fabs(x - prev_x) <= eps) {
            return x;
        }
        if (x < x2) {
            if (f_x >= f_x2) {
                x1 = x;
                f_x1 = f_x;
            } else {
                x3 = x2;
                f_x3 = f_x2;
                x2 = x;
                f_x2 = f_x;
            }
        } else {
            if (f_x >= f_x2) {
                x3 = x;
                f_x3 = f_x;
            } else {
                x1 = x2;
                f_x1 = f_x2;
                x2 = x;
                f_x2 = f_x;
            }
        }
        prev_x = x;
    }
}

/**
 * Combined Brent method of finding min value
 * @require a < b && function unimodal in [a, b]
 * @param f function for research
 * @param a left border
 * @param b right border
 * @param eps absolute accuracy
 * @return Min value in given range with given accuracy
 */
double brent(func f, double a, double c, double eps) {
    number_of_iterations = 0;
    double x, w, v, x_res, w_res, v_res, d, e;
    x = v = w = (a + c) / 2;
    x_res = w_res = v_res = f(x);
    d = e = c - a;
    while (true) {
        number_of_iterations++;
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
        if (fabs(u - x) < eps) {
            return u;
        }
        double u_res = f(u);
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

int main() {
    /**
     * absolute accuracy
     */
    double eps;

    /**
    * Function for research
    * @param x function argument
    * @return function value
    */
    auto f = [](double x) {
        return -3.0 * x * sin(0.75 * x) + exp(-2.0 * x);
    };

    std::ofstream myfile;
    myfile.open("eps.csv");
    myfile << "eps,dichotomy,golden section,fibonacci,parabolic,combined brent" << std::endl;
    for (int i = 0; i < 11; i++) {
        eps = pow(10, -i);
        std::cout << std::setprecision(i + 1);
        std::cout << std::fixed;
        myfile << eps << ",";
        dichotomy(f, 0, 2 * M_PI, eps);
        myfile << number_of_iterations << ",";
        golden_section(f, 0, 2 * M_PI, eps);
        myfile << number_of_iterations << ",";
        fibonacci(f, 0, 2 * M_PI, eps);
        myfile << number_of_iterations << ",";
        parabolic(f, 0, 2 * M_PI, eps);
        myfile << number_of_iterations << ",";
        brent(f, 0, 2 * M_PI, eps);
        myfile << number_of_iterations << std::endl;
    }
    myfile.close();

    return 0;
}