#include <fstream>
#include <iomanip>
#include "lab3.cpp"

// здесь должен быть логгер, ввод-вывод данных и запуск методов



// init csv file logger
std::ofstream logger_start(const std::string& filename, const std::string& name) {
    std::ofstream os;
    os.open(filename);
    os << std::setprecision(4) << std::fixed;
    os << name;
    return os;
}

// main log function
void log(std::ofstream& os, int n, int k, double abs_e, double rel_e) {
    os << n << "," << k << "," << abs_e << "," << rel_e << std::endl;
}

int main() {

    std::ofstream os = logger_start("test.csv", "n,k,abs,rel");
    log(os, 1, 1, 0.1, 0.1);
}