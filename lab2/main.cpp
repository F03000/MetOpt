#include<iostream>
#include <iomanip>

#include "algo.h"
#include "linear_algebra.h"

// данные функции
matrix_ A;
vector_ B;
double C;
int n;
double k;
double eps;

int number_of_iterations;

std::ofstream &operator<<(std::ofstream &s, const vector_& v) {
    s << "[";
    for (double i : v) {
        s << i << " ";
    }
    s << "]";
    return s;
}

void csv_out(std::ofstream &s, int nym, double f_x, const vector_& x, const vector_& p, const vector_& g, double a) {
    s << nym << "," << f_x << ",";
    s << x;
    s << ",";
    s << p;
    s << "," << module(p) << ",";
    s << g;
    s << ",";
    s << module(g) << "," << a << "," << eps << std::endl;
}

// Считаем квадратичную функцию в точке
double f(const vector_& x) {
    return scalar(x, (A * x) + B) + C;
}

// Градиентный спуск
double gradient_descent(const vector_& x0) {
    number_of_iterations = 0;
    vector_ x_cur = x0;
    double f_x_cur = f(x_cur);
    double alpha = k;

    while (true) {
        number_of_iterations++;
        vector_ gradient = (A * x_cur) * 2 + B;
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
double steepest_descent(const vector_& x0) {
    number_of_iterations = 0;
    vector_ x_cur = x0;
    double alpha;

    while (true) {
        number_of_iterations++;
        vector_ gradient = (A * x_cur) * 2 + B;
        if (module(gradient) < eps) {
            return f(x_cur);
        }

        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur - gradient * x);
        };

        // выбрать любой способ:
//        alpha = algo::dichotomy(f1, 0, k, eps);
//        alpha = algo::golden_section(f1, 0, k, eps);
        alpha = algo::fibonacci(f1, 0, k, eps);
//        alpha = algo::parabolic(f1, 0, k, eps);
//        alpha = algo::brent(f1, 0, k, eps);

        x_cur = x_cur - gradient * alpha;
    }
}

// Сопряженный градиент
double conjugate_gradient(const vector_& x0) {
    number_of_iterations = 0;
    vector_ x_cur = x0;
    vector_ gradient = (A * x_cur) * 2 + B;
    vector_ p_cur = vector_(n, 0) - gradient;
    double alpha;

    while (true) {
        number_of_iterations++;
        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur - gradient * x);
        };

        // выбрать любой способ:
//        alpha = algo::dichotomy(f1, 0, k, eps);
//        alpha = algo::golden_section(f1, 0, k, eps);
        alpha = algo::fibonacci(f1, 0, k, eps);
//        alpha = algo::parabolic(f1, 0, k, eps);
//        alpha = algo::brent(f1, 0, k, eps);

        vector_ x_new = x_cur + p_cur * alpha;
        vector_ gradient_new = (A * x_new) * 2 + B;
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

// Градиентный спуск
double gradient_descent_csv(const vector_& x0) {
    std::ofstream myfile;
    myfile.open("gradient_descent_n" + std::to_string(n) + "_k" + std::to_string(k) + ".csv");
    myfile << std::setprecision(4) << std::fixed;
    myfile << "n,f(x),x,p,||p||,gradient,||gradient||,alpha,eps" << std::endl;

    number_of_iterations = 0;
    vector_ x_cur = x0;
    double f_x_cur = f(x_cur);
    double alpha = k;

    while (true) {
        number_of_iterations++;
        vector_ gradient = (A * x_cur) * 2 + B;
        if (module(gradient) < eps) {
            myfile.close();
            return f_x_cur;
        }

        vector_ x_new = x_cur - (gradient * alpha);
        double f_x_new = f(x_new);

        while (f_x_new >= f_x_cur) {
            alpha /= 2;
            x_new = x_cur - (gradient * alpha);
            f_x_new = f(x_new);
        }

        vector_ p_cur = vector_ (n, 0) - gradient * alpha;
        csv_out(myfile, number_of_iterations, f(x_cur), x_cur, p_cur, gradient, alpha);

        x_cur = x_new;
        f_x_cur = f_x_new;
    }
}

// Наискорейший спуск
double steepest_descent_csv(const vector_& x0) {
    std::ofstream myfile;
    myfile.open("steepest_descent_n" + std::to_string(n) + "_k" + std::to_string(k) + ".csv");
    myfile << std::setprecision(4) << std::fixed;
    myfile << "n,f(x),x,p,||p||,gradient,||gradient||,alpha,eps" << std::endl;

    number_of_iterations = 0;
    vector_ x_cur = x0;
    double alpha;

    while (true) {
        number_of_iterations++;
        vector_ gradient = (A * x_cur) * 2 + B;
        if (module(gradient) < eps) {
            myfile.close();
            return f(x_cur);
        }

        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur - gradient * x);
        };

        // выбрать любой способ:
//        alpha = algo::dichotomy(f1, 0, k, eps);
//        alpha = algo::golden_section(f1, 0, k, eps);
        alpha = algo::fibonacci(f1, 0, k, eps);
//        alpha = algo::parabolic(f1, 0, k, eps);
//        alpha = algo::brent(f1, 0, k, eps);

        vector_ p_cur = vector_ (n, 0) - gradient * alpha;
        csv_out(myfile, number_of_iterations, f(x_cur), x_cur, p_cur, gradient, alpha);

        x_cur = x_cur - gradient * alpha;
    }
}

// Сопряженный градиент
double conjugate_gradient_csv(const vector_& x0) {
    std::ofstream myfile;
    myfile.open("conjugate_gradient_n" + std::to_string(n) + "_k" + std::to_string(k) + ".csv");
    myfile << std::setprecision(4) << std::fixed;
    myfile << "n,f(x),x,p,||p||,gradient,||gradient||,alpha,eps" << std::endl;

    number_of_iterations = 0;
    vector_ x_cur = x0;
    vector_ gradient = (A * x_cur) * 2 + B;
    vector_ p_cur = vector_(n, 0) - gradient;
    double alpha;

    while (true) {
        number_of_iterations++;
        if (module(gradient) < eps) {
            myfile.close();
            return f(x_cur);
        }
        // Одномерная оптимизация
        auto f1 = [&](double x) {
            return f(x_cur - gradient * x);
        };

        // выбрать любой способ:
//        alpha = algo::dichotomy(f1, 0, k, eps);
//        alpha = algo::golden_section(f1, 0, k, eps);
        alpha = algo::fibonacci(f1, 0, k, eps);
//        alpha = algo::parabolic(f1, 0, k, eps);
//        alpha = algo::brent(f1, 0, k, eps);

        csv_out(myfile, number_of_iterations, f(x_cur), x_cur, p_cur, gradient, alpha);

        vector_ x_new = x_cur + p_cur * alpha;
        vector_ gradient_new = (A * x_new) * 2 + B;
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
    k = 0.01;
}

// Инициализация вручную
void init() {
    n = 2;
    A = {{3, 2},
         {2, 2}};
    B = {2, -3};
    C = -3;
    k = 0.01;
}

int get_int(int k_) {
    return (std::rand() % k_);
}

// create function with specified conditions
void generate_function(int n_, int k_) {
    n = n_;

    A = matrix_(n, vector_(n, 0));
    B = vector_(n);
    for (int i = 0; i < n; ++i) {
        A[i][i] = get_int(k_);
        B[i] = get_int(k_);
    }
    C = get_int(k_);

    // число обучсловленности k_
    A[0][0] = k_;
    if (n > 1) {
        A[n - 1][n - 1] = 1;
    }
    k = 1.0 / k_;
}

void make_experiment() {
    int n_ = 11; // шаг 1
    vector_ x_0 = vector_(n, 0);
    n_ = 101; // шаг 10
    for (int k_ = 1; k_ <= n_; k_ += 10) {
        generate_function(n_, k_);
        gradient_descent_csv(x_0);
        steepest_descent_csv(x_0);
        conjugate_gradient_csv(x_0);
    }
}

// Логгер
void log(const std::string& name, double res) {
    std::cout << std::setprecision(4) << std::fixed;
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
    eps = 0.0001;
//    system("chcp 65001");

//    make_experiment();

    vector_ x_0 = vector_(n, 0);
    log("Градиентный спуск", gradient_descent(x_0));
    log("Наискорейший спуск", steepest_descent(x_0));
    log("Сопряженный градиент", conjugate_gradient(x_0));
}
