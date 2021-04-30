#include <iostream>
#include <iomanip>

#include "algo.h"
#include "linear_algebra.h"

matrix_ A;
vector_ B;
double C;
int n;
double L = 10000;
int number_of_iterations;

// Считаем квадратичную функцию в точке
double f(const vector_& x) {
    return scalar(x, (A * x) += B) + C;
}

// Градиентный спуск
double gradient_descent(const vector_& x0, double eps) {
    number_of_iterations = 0;
    vector_ x_cur = x0;
    double f_x_cur = f(x_cur);
    double alpha = L;

    while (true) {
        number_of_iterations++;
        vector_ gradient = (A * x_cur) * 2 += B;
        if (module(gradient) < eps) {
            return f_x_cur;
        }

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
double steepest_descent(const vector_& x0, double eps) {
    number_of_iterations = 0;
    vector_ x_cur = x0;
    double f_x_cur = f(x_cur);
    double alpha = L;

    while (true) {
        number_of_iterations++;
        vector_ gradient = (A * x_cur) * 2 += B;
        if (module(gradient) < eps) {
            return f_x_cur;
        }

        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur - gradient * x);
        };

        // выбрать любой способ:
//        alpha = algo::dichotomy(f1, 0, DBL_MAX, eps);
//        alpha = algo::golden_section(f1, 0, DBL_MAX, eps);
        alpha = algo::fibonacci(f1, 0, 10000, eps);
//        alpha = algo::parabolic(f1, 0, DBL_MAX, eps);
//        alpha = algo::brent(f1, 0, DBL_MAX, eps);

        x_cur -= (gradient *= alpha);
        f_x_cur = f(x_cur);
    }
}

// Сопряженный градиент
double conjugate_gradient(const vector_& x0, double eps) {
    vector_ x_cur = x0;
    vector_ gradient = (A * x_cur) * 2 += B;
    vector_ p_cur = vector_(n, 0) - gradient;
    double alpha;

    while (true) {
        number_of_iterations++;
        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur - gradient * x);
        };

        // выбрать любой способ:
//        alpha = algo::dichotomy(f1, 0, DBL_MAX, eps);
//        alpha = algo::golden_section(f1, 0, DBL_MAX, eps);
        alpha = algo::fibonacci(f1, 0, 10000, eps);
//        alpha = algo::parabolic(f1, 0, DBL_MAX, eps);
//        alpha = algo::brent(f1, 0, DBL_MAX, eps);

        vector_ x_new = x_cur + p_cur * alpha;
        vector_ gradient_new = (A * x_new) * 2 += B;
        if (module(gradient) < eps) {
            return f(x_new);
        }
        double beta;
        if (number_of_iterations % n == 0) {
            beta = 0;
        } else {
            beta = module(gradient_new) * module(gradient_new) / (module(gradient) * module(gradient));
        }
        p_cur = p_cur * beta - gradient_new;
        gradient = gradient_new;
        x_cur = x_new;
    }
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
    A = {{3, 2},
         {2, 2}};
    B = {2, -3};
    C = -3;
}

// Логгер
void log(const std::string& name, double eps, double res) {
    std::cout << "Method used:          " << name << std::endl;
    std::cout << "Absolute error        " << eps << std::endl;
    std::cout << "Result:               " << res << std::endl;
    std::cout << "Number of iterations: " << number_of_iterations << std::endl;
    std::cout << std::endl;
}

int main() {
    // Задание функции:
//    scanFunction();
    init();

    // Начальные параметры:
    double eps = 0.0001;
//    system("chcp 65001");
    std::cout << std::setprecision(5) << std::fixed;

    log("Градиентный спуск", eps, gradient_descent(vector_(n, 1), eps));
//    log("Наискорейший спуск", eps, steepest_descent(vector_(n, 1), eps));
    log("Сопряженный градиент", eps, conjugate_gradient(vector_(n, 1), eps));

}
