#ifndef METOPT_TESTS_H
#define METOPT_TESTS_H

#include "loggers.h"
#include "methods.h"

// модуль вектора
double module(const vector_& x) {
    double sum = 0;
    for (auto i : x) {
        sum += i * i;
    }
    return sqrt(sum);
}

// скалярное произведение векторов
double scalar(const vector_ &v1, const vector_ &v2) {
    double res = 0;
    for (int i = 0; i < v1.size(); ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}


// Умножение матрицы на вектор (получается вектор-столбец)
vector_ multiply(const matrix_ &m, const vector_ &v) {
    vector_ res(v.size());
    for (int i = 0; i < m.size(); i++) {
        res[i] = scalar(m[i], v);
    }
    return res;
}

// Умножение матрицы на вектор (получается вектор-столбец)
vector_ multiply(profile_matrix &m, const vector_ &v) {
    vector_ res(v.size());
    for (int k = 0; k < v.size(); k++) {
        for (int i = 0; i < v.size(); i++) {
            res[i] = 0;
            for (int j = 0; j < v.size(); j++) {
                res[i] += m.get(i, j) * v[j];
            }
        }
    }
    return res;
}

/// Тест методов на уже созданных матрицах
void simple_test() {
    matrix_ A;
    vector_ b;

    for (int i = 0; i < 5; ++i) {
        std::string input_filename = "test_" + std::to_string(i) + ".txt";
        if (!input(input_filename, A, b)) {
            continue;
        }

        std::string output_filename = "test_" + std::to_string(i) + ".csv";
        std::ofstream os = logger_start(output_filename, "");


        vector_ x;
        if (i % 2) {
            x = gauss(A, b);
        } else {
            profile_matrix P = profile_matrix(A);
            x = lu_solver(P, b);
        }

        vector_ absolute_accuracy(x.size());
        vector_ exact_solution(x.size());
        for (int j = 0; j < x.size(); j++) {
            exact_solution[j] = j + 1;
            absolute_accuracy[j] = j + 1 - x[j];
        }

        log(os, A.size(), 0, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
        os.close();
    }
}

/// Тест LU-метода на матрицах с различным числом обусловленности
void LU_test() {
    std::ofstream os = logger_start("test_LU.csv", "n,k,||x* - x_k||,||x* - x_k|| / ||x*||");

    for (int n = 10; n <= 100; n *= 10) {
        for (int k = 0; k <= 10; ++k) {
            profile_matrix p = profile_matrix(n, k);
            vector_ exact_solution(n);
            for (int i = 0; i < n; i++) {
                exact_solution[i] = i + 1;
            }
            vector_ b = multiply(p, exact_solution);
            vector_ x = lu_solver(p, b);
            vector_ absolute_accuracy(x.size());
            for (int j = 0; j < x.size(); j++) {
                absolute_accuracy[j] = exact_solution[j] - x[j];
            }
            log(os, n, k, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
        }
    }
    os.close();
}

/// Тест LU-метода на Гильбертовых матрицах
void guilbert_test() {
    std::ofstream os = logger_start("test_guilbert.csv", "n,k,||x* - x_k||,||x* - x_k|| / ||x*||");

    for (int n = 10; n <= 1000; n *= 10) {
        matrix_ g = guilbert_generator(n);
        profile_matrix p = profile_matrix(g);
        vector_ exact_solution = free_generator(n);
        vector_ b = multiply(p, exact_solution);
        vector_ x = lu_solver(p, b);
        vector_ absolute_accuracy(x.size());
        for (int j = 0; j < x.size(); j++) {
            absolute_accuracy[j] = exact_solution[j] - x[j];
        }
        log(os, n, 0, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
    }

    os.close();
}

/// Тест метода Гаусса на плотных матрицах
void gauss_test() {
    std::ofstream os = logger_start("test_gauss.csv", "n,x");
    matrix_ m;
    vector_ b;

    for (int n = 10; n <= 1000; n *= 10) {
        m = dense_generator(n);
        vector_ exact_solution = free_generator(n);
        b = multiply(m, exact_solution);
        vector_ x_gauss = gauss(m, b);
        profile_matrix pm (m);
        vector_ x_lu = lu_solver(pm, b);
        vector_ absolute_accuracy_gauss(n);
        for (int j = 0; j < n; j++) {
            absolute_accuracy_gauss[j] = exact_solution[j] - x_gauss[j];
        }
        gauss_log(os, n, x_gauss, b, m);
        vector_ absolute_accuracy_lu(n);
        for (int j = 0; j < n; j++) {
            absolute_accuracy_lu[j] = exact_solution[j] - x_lu[j];
        }
        gauss_log(os, n, x_lu, b, m);
    }
    os.close();
}

#endif //METOPT_TESTS_H
