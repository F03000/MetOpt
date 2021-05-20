#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

class profile_matrix {
    std::vector<double> di;
    std::vector<int> ia;
    std::vector<double> al;
    std::vector<double> au;

    int default_dimension = 20;
    std::vector<int> numbers_set = {0, -1, -2, -3, -4};

public:
    profile_matrix() {
        // своя матрица
        default_dimension = 3;
        di = {0, 2, 4};
        ia = {0, 0, 1, 2};
        al = {1, 3};
        au = {1, 5};
    }

    profile_matrix(int k) {
        // генерируем матрицу как в задании
        srand(time(nullptr));
        size_t size = (default_dimension * default_dimension + 1) / 2;
        al.resize(size);
        au.resize(size);
        di.resize(default_dimension);
        ia.resize(default_dimension + 1);
        int sum = 0;
        for (size_t i = 0; i < size; i++) {
            int rand = std::rand();
            al[i] = numbers_set[rand % 5];
            sum -= al[i];
            rand /= 5;
            au[i] = numbers_set[rand % 5];
            sum -= au[i];
        }
        ia[0] = 0;
        for (int i = 0; i < default_dimension; i++) {
            ia[i + 1] = ia[i] + i;
            if (i == 0) {
                di[i] = sum + pow(10, -k);
            } else {
                di[i] = sum;
            }
        }
    }

    std::vector<double> get_vector() {
        // получаем вектор f_k, как в задании
        std::vector<double> out(default_dimension);
        for (int ix = 0; ix < default_dimension; ix++)
        {
            out[ix] = 0;
            for (int jx = 0; jx < default_dimension; jx++)
                out[ix] += get(ix, jx) * (jx + 1);
        }
        return out;
    }

    int size() {
        // размер
        return default_dimension;
    }

    double get(int i, int j) {
        // получаем элемент матрицы по ее индексам
        if (i == j) {
            return di[i];
        }
        bool swapped = false;
        if (i < j) {
            std::swap(i, j);
            swapped = true;
        }
        int first = i - (ia[i + 1] - ia[i]);
        if (j < first) {
            return 0;
        } else {
            return swapped ? au[ia[i] + j - first] : al[ia[i] + j - first];
        }
    }
};
