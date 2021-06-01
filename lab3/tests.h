#ifndef METOPT_TESTS_H
#define METOPT_TESTS_H

#include "loggers.h"
#include "methods.h"
#include "linear_algebra.h"

const std::string csv_init = "n,k,||x* - x_k||,||x* - x_k|| / ||x*||\n";

/// Тест методов на уже созданных матрицах
void test_simple() {
    matrix_ A;
    vector_ b;
    vector_ exact_solution;

    for (int i = 0; i < 5; ++i) {
        std::string input_filename = "test_" + std::to_string(i) + ".txt";
        if (!input(input_filename, A, b, exact_solution)) {
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

        for (double j : x) {
            os << j << " ";
        }
        os.close();
    }
}

/// Тест LU-метода на матрицах с различным числом обусловленности
void test_lu_diagonal() {
    std::ofstream os = logger_start("test_lu_diagonal.csv", csv_init);

    for (int n = 10; n <= 100; n *= 10) {
        for (int k = 0; k <= 10; ++k) {
            profile_matrix p = profile_matrix(n, k);
            vector_ exact_solution(n);
            for (int i = 0; i < n; i++) {
                exact_solution[i] = i + 1;
            }
            vector_ b = p * exact_solution;
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
void test_lu_guilbert() {
    std::ofstream os = logger_start("test_lu_guilbert.csv", csv_init);

    for (int n = 10; n <= 1000; n *= 10) {
        matrix_ g = guilbert_generator(n);
        profile_matrix p = profile_matrix(g);
        vector_ exact_solution(n);
        for (int i = 0; i < n; i++) {
            exact_solution[i] = i + 1;
        }
        vector_ b = p * exact_solution;
        vector_ x = lu_solver(p, b);
        vector_ absolute_accuracy(x.size());
        for (int j = 0; j < x.size(); j++) {
            absolute_accuracy[j] = exact_solution[j] - x[j];
        }
        log(os, n, n * n, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
    }
    os.close();
}

/// Тест метода Гаусса на положительных матрицах с различным числом обусловленности
void test_gauss_guilbert() {
    std::ofstream os = logger_start("test_gauss_guilbert.csv", csv_init);

    for (int n = 10; n <= 1000; n *= 10) {
        matrix_ g = guilbert_generator(n);
        vector_ exact_solution(n);
        for (int i = 0; i < n; i++) {
            exact_solution[i] = i + 1;
        }
        vector_ b = g * exact_solution;
        vector_ x = gauss(g, b);
        vector_ absolute_accuracy(x.size());
        for (int j = 0; j < x.size(); j++) {
            absolute_accuracy[j] = exact_solution[j] - x[j];
        }
        log(os, n, n * n, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
    }
    os.close();
}

/// Тест метода сопряженных градиентов на уже созданных матрицах
void test_conjugate_simple() {
    matrix_ A;
    vector_ b;
    vector_ exact_solution;

    for (int i = 0; i < 5; ++i) {
        std::string input_filename = "test_" + std::to_string(i) + ".txt";
        if (!input(input_filename, A, b, exact_solution)) {
            continue;
        }

        std::string output_filename = "test_conjugate_" + std::to_string(i) + ".csv";
        std::ofstream os = logger_start(output_filename, "");

        sparse_matrix P = sparse_matrix(A);
        vector_ x = conjugate_gradient(P, b, 10e-7);

        for (double j : x) {
            os << j << " ";
        }
        os.close();
    }
}

/// Тест метода сопряженных градиентов на отрицательных матрицах с различным числом обусловленности
void test_conjugate_diagonal() {
    std::ofstream os = logger_start("test_conjugate_diagonal.csv", csv_init);

    for (int n = 10; n <= 100; n *= 10) {
        for (int k = 0; k <= 10; ++k) {
            sparse_matrix p = sparse_matrix(n, 1);
            vector_ exact_solution(n);
            for (int i = 0; i < n; i++) {
                exact_solution[i] = i + 1;
            }
            vector_ b = p * exact_solution;
            vector_ x = conjugate_gradient(p, b, 10e-7);
            vector_ absolute_accuracy(x.size());
            for (int j = 0; j < x.size(); j++) {
                absolute_accuracy[j] = exact_solution[j] - x[j];
            }
            log(os, n, k, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
        }
    }
    os.close();
}

/// Тест метода сопряженных градиентов на положительных матрицах с различным числом обусловленности
void test_conjugate_reverse_diagonal() {
    std::ofstream os = logger_start("test_conjugate_reverse_diagonal.csv", csv_init);

    for (int n = 10; n <= 100; n *= 10) {
        for (int k = 0; k <= 10; ++k) {
            sparse_matrix p = sparse_matrix(n, -1);
            vector_ exact_solution(n);
            for (int i = 0; i < n; i++) {
                exact_solution[i] = i + 1;
            }
            vector_ b = p * exact_solution;
            vector_ x = conjugate_gradient(p, b, 10e-7);
            vector_ absolute_accuracy(x.size());
            for (int j = 0; j < x.size(); j++) {
                absolute_accuracy[j] = exact_solution[j] - x[j];
            }
            log(os, n, k, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
        }
    }
    os.close();
}

/// Тест метода сопряженных градиентов на Гильбертовых матрицах
void test_conjugate_guilbert() {
    std::ofstream os = logger_start("test_conjugate_guilbert.csv", csv_init);

    for (int n = 10; n <= 1000; n *= 10) {
        matrix_ g = guilbert_generator(n);
        sparse_matrix p = sparse_matrix(g);
        vector_ exact_solution(n);
        for (int i = 0; i < n; i++) {
            exact_solution[i] = i + 1;
        }
        vector_ b = g * exact_solution;
        vector_ x = conjugate_gradient(p, b, 10e-7);
        vector_ absolute_accuracy(x.size());
        for (int j = 0; j < x.size(); j++) {
            absolute_accuracy[j] = exact_solution[j] - x[j];
        }
        log(os, n, n * n, module(absolute_accuracy), module(absolute_accuracy) / module(exact_solution));
    }
    os.close();
}

#endif //METOPT_TESTS_H
