#ifndef METOPT_LINEAR_ALGEBRA_H
#define METOPT_LINEAR_ALGEBRA_H

#include <cmath>
#include <vector>

typedef std::vector<std::vector<double>> matrix_;
typedef std::vector<double> vector_;

// скалярное произведение векторов
double scalar(const vector_& v1, const vector_& v2) {
    double res = 0;
    for (int i = 0; i < v1.size(); ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

// Умножение вектора на число
vector_ operator*=(vector_ &v, double a) {
    for (double & i : v) {
        i *= a;
    }
    return v;
}

vector_ operator*(const vector_ &v, double a) {
    vector_ t = vector_(v.size());
    for (int i = 0; i < v.size(); ++i) {
        t[i] = v[i] * a;
    }
    return t;
}

// Сложение векторов
vector_ operator+=(vector_ v1, const vector_& v2) {
    for (int i = 0; i < v1.size(); i++) {
        v1[i] += v2[i];
    }
    return v1;
}

vector_ operator+(const vector_& v1, const vector_& v2) {
    return vector_(v1.size(), 0) += v1 += v2;
}

// Вычитание векторов
vector_ operator-=(vector_ v1, const vector_& v2) {
    for (int i = 0; i < v1.size(); i++) {
        v1[i] -= v2[i];
    }
    return v1;
}

vector_ operator-(const vector_& v1, const vector_& v2) {
    return vector_(v1.size(), 0) += v1 -= v2;
}

// Умножение вектора на матрицу (получается вектор-строка)
vector_ operator*(const vector_& v, const matrix_& m) {
    vector_ result(v.size(), 0);
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            result[i] += m[j][i] * v[j];
        }
    }
    return result;
}

// Умножение матрицы на вектор (получается вектор-столбец)
vector_ operator*(const matrix_& m, const vector_& v) {
    vector_ res(v.size());
    for (int i = 0; i < m.size(); i++) {
        res[i] = scalar(m[i], v);
    }
    return res;
}

// Умножение матрицы на матрицу
matrix_ operator*(const matrix_& m1, const matrix_& m2) {
    matrix_ res(m1.size());
    for (int i = 0; i < res.size(); i++) {
        res[i] = m1[i] * m2;
    }
    return res;
}

// Умножение матрицы на число
matrix_ operator*=(matrix_ m, double a) {
    for (vector_ &v: m) {
        v = v * a;
    }
    return m;
}

matrix_ operator*(matrix_ m, double a) {
    matrix_ t = matrix_(m.size(), vector_(m.size()));
    for (int i = 0; i < t.size(); ++i) {
        for (int j = 0; j < t.size(); ++j) {
            t[i][j] = m[i][j] * a;
        }
    }
    return t;
}

// длина вектора
double module(const vector_& v) {
    return sqrt(scalar(v, v));
}

// Дистанция между точками
double euclidean_distance(const vector_& v1, const vector_& v2) {
    return module(v1 - v2);
}

#endif //METOPT_LINEAR_ALGEBRA_H
