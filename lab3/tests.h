#ifndef METOPT_TESTS_H
#define METOPT_TESTS_H

#include "loggers.h"
#include "methods.h"

/// Тест методов на уже созданных матрицах
void simple_test() {  // TODO: написать больше тестов
    matrix_ A;
    vector_ b;

    for (int i = 0; i < 1; ++i) {
        std::string root = "/home/rytuo/work/metOpt/MetOpt/lab3/";
        std::string input_filename = root + "tests/test" + std::to_string(i) + ".txt";
        std::string output_filename = root + "output/test" + std::to_string(i) + ".csv";

        if (!input(input_filename, A, b)) {
            continue;
        }

        profile_matrix P = profile_matrix(A);
        vector_ x = lu_solving(P, b);

//        vector_ x = gauss(A, b);

        std::ofstream os;
        os.open(output_filename);
        for (double j : x) {
            os << j << " ";
        }
        os.close();
    }
}

/// Тест LU-метода на матрицах с различным числом обусловленности
void LU_test() {

}

/// Тест LU-метода на Гильбертовых матрицах
void guilbert_test() {

}

/// Тест метода Гаусса на плотных матрицах
void gauss_test() {

}

#endif //METOPT_TESTS_H
