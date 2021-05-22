#ifndef METOPT_TESTS_H
#define METOPT_TESTS_H

#include "loggers.h"
#include "methods.h"

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

        // TODO: вычислить погрешность ????
        log(os, A.size(), 0, 0.1, 0.1);
        os.close();
    }
}

/// Тест LU-метода на матрицах с различным числом обусловленности
void LU_test() {
    std::ofstream os = logger_start("test_LU.csv","n,k,||x* - x_k||,||x* - x_k|| / ||x*||");

    for (int n = 10; n <= 1000; n *= 10) {
        for (int k = 0; k <= 10; ++k) {
            profile_matrix p = profile_matrix(n, k);
            vector_ b = free_generator(n);
            vector_ x = lu_solver(p, b);

            // TODO: вычислить погрешность
            log(os, n, k, 0.1, 0.1);
        }
    }

    os.close();
}

/// Тест LU-метода на Гильбертовых матрицах
void guilbert_test() {
    std::ofstream os = logger_start("test_guilbert.csv","n,k,||x* - x_k||,||x* - x_k|| / ||x*||");

    for (int n = 10; n <= 1000; n *= 10) {
        matrix_ g = guilbert_generator(n);
        profile_matrix p = profile_matrix(g);
        vector_ b = free_generator(n);
        vector_ x = lu_solver(p, b);

        // TODO: вычислить погрешность
        log(os, n, 0, 0.1, 0.1);
    }

    os.close();
}

/// Тест метода Гаусса на плотных матрицах
void gauss_test() { // TODO: как-то сравнить с LU
    std::ofstream os = logger_start("test_gauss.csv", "n,x");
    matrix_ m;
    vector_ b;

    for (int n = 10; n <= 1000; n *= 10) {
        m = dense_generator(n);
        b = free_generator(n);
        vector_ x = gauss(m, b);
        gauss_log(os, n, x, b, m);
    }

    os.close();
}

#endif //METOPT_TESTS_H
