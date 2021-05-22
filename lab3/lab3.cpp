#include "lab3.h"

// здесь должен быть код методов и различные эксперименты

/**
 * Алгоритм решения слау на основе LU-разложения
 * На вход подается квадратная матрица ?профильного? формата
 * На выходе вектор x - одно из решений слау (если есть)
 * */
std::vector<double> lu_solving(profile_matrix &A, vector_ &b) {
    int n = A.size();
    for (int j = 0; j < n; j++) {
        for (int i = j + 1; i < n; i++) {
            A.set(i, j, A.get(i, j) / A.get(j, j));
            for (int k = j + 1; k < n; k++) {
                A.set(i, k, A.get(i, k) - A.get(i, j) * A.get(j, k));
            }
        }
    }
    vector_ y = vector_(n);
    // решить Ly = b прямым ходом
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= A.get(i, j) * y[j];
        }
    }

    vector_ x = vector_(n);
    // решить Ux = y обратным ходом
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A.get(i, j) * x[j];
        }
        x[i] /= A.get(i, i);
    }
    return x;
}

/**
 * Алгоритм решеия слау двупроходным методом Гаусса с выбором ведущего элемента
 * На вход подается квадратная матрица ?плотного? формата
 * На выходе вектор x - одно из решений слау (если есть)
 * */
std::vector<double> gauss(matrix_ &A, vector_ &b) {
    int n = (int)A.size();

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

// TODO: написать генератор матриц с разными обусловленностями по схеме из условия лабы

// TODO: написать генератор гильбертовых матриц по схеме из условия лабы