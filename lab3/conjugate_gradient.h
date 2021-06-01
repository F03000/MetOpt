#ifndef METOPT_CONJUGATE_GRADIENT_H
#define METOPT_CONJUGATE_GRADIENT_H

#include "linear_algebra.h"

// данные функции
int number_of_iterations;

// Считаем квадратичную функцию в точке
double f(matrix_ &A, vector_ &B, double C, const vector_& x) {
    return scalar(x, (A * x) + B) + C;
}

// Сопряженный градиент
vector_ conjugate_gradient(sparse_matrix &A, vector_ &B, double C, const vector_ &x0, double eps) {
    int n = (int)A.size();
    number_of_iterations = 0;
    vector_ x_cur = x0;
    vector_ gradient = (A * x_cur) * 2 + B;
    vector_ p_cur = vector_(n, 0) - gradient;
    double alpha;

    while (true) {
        if (module(gradient) < eps) {
            return x_cur;
        }
        number_of_iterations++;

        // Одномерная оптимизация
        alpha = -((scalar(gradient, p_cur)) / (scalar(A * p_cur, p_cur))) / 2;

//        auto f1 = [&](double x) {
//            return f(A, B, C, x_cur + p_cur * x);
//        };
//        alpha = algo::fibonacci(f1, 0, max_alpha, eps);

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

#endif //METOPT_CONJUGATE_GRADIENT_H
