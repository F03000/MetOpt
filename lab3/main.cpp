#include <fstream>
#include <iomanip>
#include "lab3.cpp"

// здесь должен быть логгер, ввод-вывод данных и запуск методов

void input(const std::string &filename, matrix_ &A, vector_ &b) {
    std::ifstream is;
    is.open(filename);
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
    input("input.txt", A, b);

    std::ofstream os = logger_start("output.csv", "n,k,abs,rel");
    log(os, 1, 1, 0.1, 0.1);
    os.close()
}