#ifndef METOPT_LOGGERS_H
#define METOPT_LOGGERS_H

#include "matrix.h"

#include <fstream>
#include <iostream>
#include <iomanip>

std::string ROOT = "/home/rytuo/work/metOpt/MetOpt/lab3/";

/// Ввод из файла: размерность n, матрица A, вектор b
int input(const std::string &filename, matrix_ &A, vector_ &b) {
    std::ifstream is(ROOT + "test/" + filename);
    if (!is) {
        std::cerr << "Could not open file: " + filename << std::endl;
        return 0;
    } else {
        std::cout << "New test: " + filename << std::endl;
    }
    int n;
    is >> n;
    A = matrix_(n, vector_(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            is >> A[i][j];
        }
    }
    b = vector_(n);
    for (int i = 0; i < n; ++i) {
        is >> b[i];
    }
    is.close();
    return 1;
}

/// Инициализация логгера в выбранный файл
std::ofstream logger_start(const std::string &filename, const std::string &name) {
    std::ofstream os;
    os.open(ROOT + "output/" + filename);
    os << std::setprecision(4) << std::fixed;
    os << name;
    return os;
}

/// Логирующая функция
void log(std::ofstream &os, int n, int k, double abs_e, double rel_e) {
    os << n << "," << k << "," << abs_e << "," << rel_e << std::endl;
}

/// Логирующая функция для метода Гаусса
void gauss_log(std::ofstream &os, int n, const vector_ &x, const vector_ &b, const matrix_ &m) {
    os << n << ",[";
    for (double i : x) {
        os << i << " ";
    }
    os << "]" << std::endl;
}

#endif //METOPT_LOGGERS_H
