#include <iostream>
#include <utility>

#include "linear_algebra.h"
#include "algo.h"

matrix A;
vector B;
double C;
int n;

// Считаем квадратичную функцию в точке
double f(const vector& x) {
    return scalar(x, (A * x) += B) + C;
}

// Градиентный спуск
double gradient_descent(vector x0, double alpha, double eps) {
    vector x_cur = std::move(x0);
    double f_x_cur = f(x_cur);

    while (true) {
        vector gradient = (A * x_cur) += B;
        if (module(gradient) < eps) {
            return f_x_cur;
        }

        vector x_new = x_cur - (gradient * alpha);
        double f_x_new = f(x_new);
        while (f_x_new >= f_x_cur) {
            alpha /= 2;
            x_new = x_cur - (gradient * alpha);
            f_x_new = f(x_new);
        }
        x_cur = x_new;
        f_x_cur = f_x_new;
    }
}

// Наискорейший спуск
double steepest_descent(const vector& x0, double alpha, double eps) {
    const vector& x_cur = x0;
    double f_x_cur = f(x_cur);

    while (true) {
        vector gradient = (A * x_cur) += B;
        if (module(gradient) < eps) {
            return f_x_cur;
        }

        // single-dimension minimisation
        auto f1 = [&](double x) {
            return f(x_cur * x);
        };
//        alpha = algo::dichotomy(f1, 0, alpha, eps);
//        alpha = algo::golden_section(f1, 0, alpha, eps);
        alpha = algo::fibonacci(f1, 0, alpha, eps);
//        alpha = algo::parabolic(f1, 0, alpha, eps);
//        alpha = algo::brent(f1, 0, alpha, eps);

        x_cur -= (gradient *= alpha);
        f_x_cur = f(x_cur);
    }
}

// Сопряженный градиент
double conjugate_gradient(vector x0, double alpha, double eps) {
    // TODO
    return 0;
}

// Ввод функции, функция задается в виде:
// - матрицы A (матричная форма квадратичной функции),
// - вектора B (соответствующие коэффициенты при x),
// - числа С (свободного члена)
void scanFunction() {
    std::cout << "Enter dimension n: " << std::endl;
    std::cin >> n;
    B.resize(n);
    A.resize(n, vector(n));
    std::cout << "Enter matrix A" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
        }
    }
    std::cout << "Enter vector B" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cin >> B[i];
    }
    std::cout << "Enter number C" << std::endl;
    std::cin >> C;
}

void init() {
    n = 2;
    A = {{6, 2},
         {2, 4}};
    B = {2, -3};
    C = -3;
}

int main() {
    // choose initialisation:
//    scanFunction();
    init();

    std::cout << gradient_descent(vector(n, 1), 0.1, 0.01);
    std::cout << steepest_descent(vector(n, 1), 0.1, 0.01);
    std::cout << conjugate_gradient(vector(n, 1), 0.1, 0.01);
}
