#include <iostream>
#include <cmath>

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


// TODO: fibonacci method
double fibonacci() {

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
    if (std::abs(x - prev_x) <= eps) {
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

// TODO: combined Brent method
double brent() {

}

/**
 * Console logger function
 * @param name: name of method
 * @param eps: accuracy
 * @param res: result
 */
void log(const std::string& name, double eps, double res) {
    std::cout << "Method used:          " << name << std::endl;
    std::cout << "Absolute error        " << eps << std::endl;
    std::cout << "Result:               " << res << std::endl;
    std::cout << "Number of iterations: " << number_of_iterations << std::endl;
    std::cout << std::endl;
}

int main() {
    double eps = 10e-5;
    log("Dichotomy", eps, dichotomy(0, M_2_PI, eps));
    log("Golden section", eps, golden_section(0, M_2_PI, eps));
//    log("Fibonacci", eps, );
//    log("Parabolic", eps, );
//    log("Combined Brent", eps, );
    return 0;
}