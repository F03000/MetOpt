#include <iostream>
#include <utility>

#include "algo.h"
#include "linear_algebra.h"

matrix_ A;
vector_ B;
double C;
int n;
int number_of_iterations;

// Считаем квадратичную функцию в точке
double f(const vector_& x) {
    return scalar(x, (A * x) += B) + C;
}

// Градиентный спуск
double gradient_descent(vector_ x0, double alpha, double eps) {
    number_of_iterations = 0;
    vector_ x_cur = std::move(x0);
    double f_x_cur = f(x_cur);

    while (true) {
        vector_ gradient = (A * x_cur) += B;
        if (module(gradient) < eps) {
            return f_x_cur;
        }

        number_of_iterations++;
        vector_ x_new = x_cur - (gradient * alpha);
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
double steepest_descent(const vector_& x0, double alpha, double eps) {
    number_of_iterations = 0;
    const vector_& x_cur = x0;
    double f_x_cur = f(x_cur);

    while (true) {
        vector_ gradient = (A * x_cur) += B;
        if (module(gradient) < eps) {
            return f_x_cur;
        }

        number_of_iterations++;
        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur * x);
        };

        // выбрать любой способ:
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
double conjugate_gradient(vector_ x0, double alpha, double eps) {
    // TODO
    return 0;
}

// Консольный ввод функции, задается в виде:
// - матрицы A (матричная форма квадратичной функции),
// - вектора B (соответствующие коэффициенты при x),
// - числа С (свободного члена)
void scanFunction() {
    std::cout << "Enter dimension n: " << std::endl;
    std::cin >> n;
    B.resize(n);
    A.resize(n, vector_(n));
    std::cout << "Enter matrix_ A" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A[i][j];
        }
    }
    std::cout << "Enter vector_ B" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cin >> B[i];
    }
    std::cout << "Enter number C" << std::endl;
    std::cin >> C;
}

// Инициализация вручную
void init() {
    n = 2;
    A = {{211, -420},
         {0, 211}};
    B = {-192, 50};
    C = -25;
}

// Логгер
void log(const std::string& name, double eps, double res) {
    std::cout << "Method used:          " << name << std::endl;
    std::cout << "Absolute error        " << eps << std::endl;
    std::cout << "Result:               " << res << std::endl;
    std::cout << "Number of iterations: " << algo::number_of_iterations << std::endl;
    std::cout << std::endl;
}

int main() {
    // Задание функции:
//    scanFunction();
    init();

    // Начальные параметры:
    double eps = 0.001;
    double alpha = 0.1;

    log("Градиентный спуск", eps, gradient_descent(vector_(n, 1), alpha, eps));
    log("Наискорейший спуск", eps, steepest_descent(vector_(n, 1), alpha, eps));
//    log("Сопряженный градиент", eps, conjugate_gradient(vector_(n, 1), alpha, eps));
}
