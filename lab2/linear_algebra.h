#ifndef METOPT_LINEAR_ALGEBRA_H
#define METOPT_LINEAR_ALGEBRA_H

#include <cmath>

typedef std::vector<std::vector<double>> matrix;
typedef std::vector<double> vector;

// скалярное произведение векторов
double scalar(const vector& v1, const vector& v2) {
    double res = 0;
    for (int i = 0; i < v1.size(); ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

// Умножение вектора на число
vector operator*=(vector &v, double a) {
    for (double & i : v) {
        i *= a;
    }
    return v;
}

vector operator*(const vector &v, double a) {
    vector t = vector(v.size());
    for (int i = 0; i < v.size(); ++i) {
        t[i] = v[i] * a;
    }
    return t;
}

// Сложение векторов
vector operator+=(vector v1, const vector& v2) {
    for (int i = 0; i < v1.size(); i++) {
        v1[i] += v2[i];
    }
    return v1;
}

vector operator+(const vector& v1, const vector& v2) {
    return vector(v1.size(), 0) += v1 += v2;
}

// Вычитание векторов
vector operator-=(vector v1, const vector& v2) {
    for (int i = 0; i < v1.size(); i++) {
        v1[i] -= v2[i];
    }
    return v1;
}

vector operator-(const vector& v1, const vector& v2) {
    return vector(v1.size(), 0) += v1 -= v2;
}

// Умножение вектора на матрицу (получается вектор-строка)
vector operator*(const vector& v, const matrix& m) {
    vector result(v.size(), 0);
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            result[i] += m[j][i] * v[j];
        }
    }
    return result;
}

// Умножение матрицы на вектор (получается вектор-столбец)
vector operator*(const matrix& m, const vector& v) {
    vector res(v.size());
    for (int i = 0; i < m.size(); i++) {
        res[i] = scalar(m[i], v);
    }
    return res;
}

// Умножение матрицы на матрицу
matrix operator*(const matrix& m1, const matrix& m2) {
    matrix res(m1.size());
    for (int i = 0; i < res.size(); i++) {
        res[i] = m1[i] * m2;
    }
    return res;
}

// Умножение матрицы на число
matrix operator*=(matrix m, double a) {
    for (vector &v: m) {
        v = v * a;
    }
    return m;
}

matrix operator*(matrix m, double a) {
    matrix t = matrix(m.size(), vector(m.size()));
    for (int i = 0; i < t.size(); ++i) {
        for (int j = 0; j < t.size(); ++j) {
            t[i][j] = m[i][j] * a;
        }
    }
    return t;
}

// длина вектора
double module(const vector& v) {
    return sqrt(scalar(v, v));
}

// Дистанция между точками
double euclidean_distance(const vector& v1, const vector& v2) {
    return module(v1 - v2);
}

#endif //METOPT_LINEAR_ALGEBRA_H
