#include <iostream>
#include <vector>
#include <cmath>

std::vector<std::vector<double>> A_matrix;
std::vector<double> B_vector;
double c;
int n;

// Умножение матрицы на вектор (получается столбец)
std::vector<double> multiply(std::vector<std::vector<double>> matrix, std::vector<double> vector) {
    std::vector<double> result(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

// Умножение вектора на матрицу (получается строка)
std::vector<double> multiply(std::vector<double> vector, std::vector<std::vector<double>> matrix) {
    std::vector<double> result(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[j][i] * vector[j];
        }
    }
    return result;
}

// Умножение матрицы на матрицу
std::vector<std::vector<double>> multiply(
        std::vector<std::vector<double>> matrix_first,
        std::vector<std::vector<double>> matrix_second) {
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += matrix_first[i][k] * matrix_second[k][j];
            }
        }
    }
    return result;
}

// Умножение матрицы на число
std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> matrix, double number) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] *= number;
        }
    }
    return matrix;
}

// Умножение вектора на число
std::vector<double> multiply(std::vector<double> vector, double number) {
    for (int i = 0; i < n; i++) {
        vector[i] *= number;
    }
    return vector;
}

// Вычитание векторов
std::vector<double> subtract(std::vector<double> vector_first, std::vector<double> vector_second) {
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = vector_first[i] - vector_second[i];
    }
    return result;
}

// Сложение векторов
std::vector<double> add(std::vector<double> vector_first, std::vector<double> vector_second) {
    std::vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = vector_first[i] + vector_second[i];
    }
    return result;
}

double module(std::vector<double> vector) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += vector[i] * vector[i];
    }
    return sqrt(result);
}

// Дистанция между точками
double euclidean_distance(std::vector<double> vector_first, std::vector<double> vector_second) {
    double tmp = 0;
    for (int i = 0; i < n; i++) {
        tmp += (vector_first[i] - vector_second[i]) * (vector_first[i] - vector_second[i]);
    }
    return sqrt(tmp);
}

// Ввод функции, функция задается в виде матрицы A (гугли матричная форма квадратичной функции),
// вектора B - B[i] коэффициент при x[i] и числа С - свободного члена
void scanFunction() {
    std::cout << "Enter dimension n: " << std::endl;
    std::cin >> n;
    B_vector.resize(n);
    A_matrix.resize(n, std::vector<double>(n));
    std::cout << "Enter matrix A" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> A_matrix[i][j];
        }
    }
    std::cout << "Enter vector B" << std::endl;
    for (int i = 0; i < n; i++) {
        std::cin >> B_vector[i];
    }
    std::cout << "Enter number C" << std::endl;
    std::cin >> c;
}

// Считаем функцию в точке
double function(std::vector<double> x) {
    std::vector<std::vector<double>> tmp_matrix(n, std::vector<double>(n, 0));
    std::vector<double> result_vector;
    for (int i = 0; i < n; i++) {
        tmp_matrix[i][i] = x[i];
    }
    tmp_matrix = multiply(A_matrix, tmp_matrix);
    tmp_matrix = multiply(tmp_matrix, 0.5);
    result_vector = multiply(x, tmp_matrix);
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += result_vector[i];
        result += B_vector[i] * x[i];
    }
    result += c;
    return result;
}

// Градиентный спуск
double gradient_descent(std::vector<double> beginning_point, double learning_rate, double eps) {
    std::vector<double> x_current = std::move(beginning_point), x_previous;
    double x_value_current = function(x_current);
    while (true) {
        std::vector<double> gradient = add(multiply(A_matrix, x_current), B_vector);
        if (module(gradient) < eps) {
            return x_value_current;
        }
        std::vector<double> x_new = subtract(x_current, multiply(gradient, learning_rate));
        double x_value_new = function(x_new);
        if (x_value_new < x_value_current) {
            x_current = x_new;
            x_value_current = x_value_new;
        } else {
            learning_rate /= 2;
        }
    }
}

int main() {
//    scanFunction();
    n = 2;
    A_matrix = {{6, 2},
                {2, 4}};
    B_vector = {2, -3};
    c = -3;

    std::cout << gradient_descent(std::vector<double>(n, 1), 0.1, 0.01);


}
