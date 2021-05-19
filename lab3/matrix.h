#ifndef METOPT_MATRIX_H
#define METOPT_MATRIX_H

#include <utility>
#include <vector>

class profile_matrix {
private:
    std::vector<double> d;
    std::vector<int> ia;
    std::vector<double> al;
    std::vector<double> au;
public:

};

class default_matrix {
private:
    std::vector<std::vector<double>> m;
public:
    explicit default_matrix(std::vector<std::vector<double>> m) {
        this->m = std::move(m);
    }
};


#endif //METOPT_MATRIX_H
