#include <iostream>
#include <vector>


void mul_vec_num(std::vector<double>& v, double num) {
    size_t N = v.size();
    for (size_t i = 0; i < N; i++) {
        v[i] = v[i] * num;
    }
}


std::vector<double> sum_two_vec(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> v;
    size_t N = v1.size();
    
    if (v1.size() != v2.size()) {
        throw std::runtime_error("Разные размеры векторов");
    } else {
        for (size_t i = 0; i < N; i++) {
            v[i] = v1[i] + v2[i];
        }
    }
    return v;
}


double mul_scal(std::vector<double>& v1, std::vector<double>& v2) {
    double res = 0.0;
    size_t N = v1.size();
    
    for (size_t i = 0; i < N; i++) {
        res += v1[i] * v2[i];
    }

    return res;
}