#ifndef METOPT_LAB3_H
#define METOPT_LAB3_H

// здесь должны быть инклуды, дефайны и вспомогательные структуры для методов

#include <utility>
#include <vector>

typedef std::vector<double> vector_;
typedef std::vector<std::vector<double>> matrix_;

class profile_matrix {
private:
    std::vector<double> d;
    std::vector<int> ia;
    std::vector<double> al;
    std::vector<double> au;
public:
    /// Передача массивов по ссылке
    explicit profile_matrix(vector_ &d, std::vector<int> &ia, vector_ &al, vector_ &au) {
        this->d = d;
        this->ia = ia;
        this->al = al;
        this->au = al;
    }
};

#endif //METOPT_LAB3_H
