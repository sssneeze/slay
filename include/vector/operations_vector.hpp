#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Vec {
    void mul_num(std::vector<double>& v, double num);
    std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> sub(const std::vector<double>& a, const std::vector<double>& b);
    double dot(const std::vector<double>& a, const std::vector<double>& b);
    double norm(const std::vector<double>& v);
    double norm_inf(const std::vector<double>& v);
}