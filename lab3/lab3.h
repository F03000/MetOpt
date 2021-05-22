#ifndef METOPT_LAB3_H
#define METOPT_LAB3_H

// здесь должны быть инклуды, дефайны и вспомогательные структуры для методов

#include <utility>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

typedef std::vector<double> vector_;
typedef std::vector<std::vector<double>> matrix_;

class profile_matrix {
private:
    int n;
    vector_ d, al, au;
    std::vector<int> il, iu;

public:
    explicit profile_matrix(matrix_ &m) {
        n = (int)m.size();

        d = vector_ ();
        al = vector_();
        au = vector_();
        il = std::vector<int>();
        iu = std::vector<int>();
        for (int i = 0, j; i < n; ++i) {
            d.push_back(m[i][i]);

            il.push_back((int)al.size());
            for (j = 0; j < i && m[i][j] == 0; ++j);
            for (; j < i; ++j) { al.push_back(m[i][j]); }

            iu.push_back((int)au.size());
            for (j = 0; j < i && m[j][i] == 0; ++j);
            for (; j < i; ++j) { au.push_back(m[j][i]); }
        }
        il.push_back((int)al.size());
        iu.push_back((int)au.size());
    }

    int size() const {
        return n;
    }

    /// получаем элемент матрицы по индексу
    double get(int i, int j) {
        if (i == j) {
            return d[i];
        } else if (i > j) {
            // номер относительно начала массива
            int k = j - i + il[i + 1] - il[i];
            return k < 0 ? 0 : al[k];
        } else {
            // номер относительно начала массива
            int k = i - j + iu[j + 1] - iu[j];
            return k < 0 ? 0 : au[k];
        }
    }

    /// меняем элемент по индексу
    void set(int i, int j, double v) {
        if (i == j) {
            d[i] = v;
        } else if (i > j) {
            // номер относительно начала массива
            int k = j - i + il[i + 1] - il[i];
            if (k < 0) {
                for (int l = 0; l < -k - 1; ++l) {
                    al.insert(al.begin() + il[i], 0);
                }
                al.insert(al.begin() + il[i], v);
                for (int l = i + 1; l < il.size(); ++l) {
                    il[l] -= k;
                }
            } else {
                al[k] = v;
            }
        } else {
            // номер относительно начала массива
            int k = i - j + iu[j + 1] - iu[j];
            if (k < 0) {
                for (int l = 0; l < -k - 1; ++l) {
                    au.insert(au.begin() + iu[j], 0);
                }
                au.insert(au.begin() + iu[j], v);
                for (int l = j + 1; l < iu.size(); ++l) {
                    iu[l] -= k;
                }
            } else {
                au[k] = v;
            }
        }
    }
};


#endif //METOPT_LAB3_H
