#include "lab3.h"

// здесь должен быть код методов и различные эксперименты

/**
 * Алгоритм решения слау на основе LU-разложения
 * На вход подается квадратная матрица ?профильного? формата
 * На выходе вектор x - одно из решений слау (если есть)
 * */
std::vector<double> lu_solving(matrix_ &A, vector_ &b) {
    int n = (int)A.size();

    // LU-разложение матрицы A
    // TODO: хранить новые матрицы на месте предыдущей
    matrix_ L = matrix_(n, vector_(n, 0));
    matrix_ U = matrix_(n, vector_(n, 0));
    L[0][0] = A[0][0];
    U[0][0] = 1;
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            L[i][j] = A[i][j];
            U[j][i] = A[j][i];
            for (int k = 0; k < j; ++k) {
                L[i][j] -= L[i][k] * U[k][j];
                U[j][i] -= L[j][k] * U[k][i];
            }
            U[j][i] /= L[j][j];
        }

        L[i][i] = A[i][i];
        for (int k = 0; k < i; ++k) {
            L[i][i] -= L[i][k] * U[k][i];
        }

        U[i][i] = 1;
    }

    // решить Ly = b прямым ходом
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            b[j] -= b[i] * L[j][i] / L[i][i];
        }
    }
    vector_ y = vector_(n);
    for (int i = 0; i < n; ++i) {
        y[i] = b[i] / L[i][i];
    }

    // решить Ux = y обратным ходом
    for (int i = n - 1; i > 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            y[j] -= y[i] * U[j][i];
        }
    }
    return y;
}

/**
 * Алгоритм решеия слау двупроходным методом Гаусса с выбором ведущего элемента
 * На вход подается квадратная матрица ?плотного? формата
 * На выходе вектор x - одно из решений слау (если есть)
 * */
std::vector<double> gauss(matrix_ &A, vector_ &b) {
    int n = (int)A.size();
    // TODO: матрицу обрабатывать на месте, желательно умно менять строчки местами (менять ссылки???)

    // прямой ход
    for (int i = 0; i < n - 1; ++i) {
        // выбор опорного элемента
        double m = std::abs(A[i][i]);
        int ind = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(A[j][i]) > m) {
                m = std::abs(A[j][i]);
                ind = j;
            }
        }
        if (m < 10e-7) {
            return vector_(n, 0);
        }
        std::swap(A[i], A[ind]);
        std::swap(b[i], b[ind]);

        // уменьшаем все строчки
        for (int j = i + 1; j < n; ++j) {
            double q = A[j][i] / A[i][i];
            A[j][i] = 0;
            for (int k = i + 1; k < n; ++k) {
                A[j][k] -= q * A[i][k];
            }
            b[j] -= q * b[i];
        }
    }

    // обратный ход
    for (int i = n - 1; i > 0; --i) {
        // уменьшаем все строчки
        for (int j = i - 1; j >= 0; --j) {
            double q = A[j][i] / A[i][i];
            A[j][i] = 0;
            for (int k = i + 1; k < n; ++k) {
                A[j][k] -= q * A[i][k];
            }
            b[j] -= q * b[i];
        }
    }

    vector_ x = vector_(n);
    for (int i = 0; i < n; ++i) {
        x[i] = b[i] / A[i][i];
    }

    return x;
}