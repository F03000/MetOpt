#include <fstream>
#include <iostream>
#include <iomanip>
#include "lab3.cpp"

// здесь должен быть логгер, ввод-вывод данных и запуск методов

/// Ввод из файла: размерность n, матрица A, вектор b
int input(const std::string &filename, matrix_ &A, vector_ &b) {
    std::ifstream is(filename);
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
    os.open(filename);
    os << std::setprecision(4) << std::fixed;
    os << name;
    return os;
}

/// Логирующая функция
void log(std::ofstream& os, int n, int k, double abs_e, double rel_e) {
    os << n << "," << k << "," << abs_e << "," << rel_e << std::endl;
}

int main() {
    matrix_ A;
    vector_ b;

    // TODO: написать больше тестов
    for (int i = 0; i < 1; ++i) {
        std::string root = "/home/rytuo/work/metOpt/MetOpt/lab3/";
        std::string input_filename = root + "tests/test" + std::to_string(i) + ".txt";
        std::string output_filename = root + "output/test" + std::to_string(i) + ".csv";

        if (!input(input_filename, A, b)) {
            continue;
        }

        // TODO: Передавать логгер функцию в методы
        // TODO: Тщательно задебажить функции
//        vector_ x = lu_solving(A, b);
        vector_ x = gauss(A, b);

        std::ofstream os;
        os.open(output_filename);
        for (double j : x) {
            os << j << " ";
        }
        os.close();
    }
}