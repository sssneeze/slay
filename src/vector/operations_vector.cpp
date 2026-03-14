#include <iostream>
#include <cmath>
#include <stdexcept>
#include "vector/operations_vector.hpp"

namespace Vec {
    void mul_num(std::vector<double>& v, double num) {
        for (auto& val : v) val *= num;
    }

    std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) 
            throw std::runtime_error("Vec::add: size mismatch");
        std::vector<double> res(a.size());
        for (size_t i = 0; i < a.size(); ++i) 
            res[i] = a[i] + b[i];
        return res;
    }

    std::vector<double> sub(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) 
            throw std::runtime_error("Vec::sub: size mismatch");
        std::vector<double> res(a.size());
        for (size_t i = 0; i < a.size(); ++i) 
            res[i] = a[i] - b[i];
        return res;
    }

    double dot(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) 
            throw std::runtime_error("Vec::dot: size mismatch");
        double res = 0.0;
        for (size_t i = 0; i < a.size(); ++i) 
            res += a[i] * b[i];
        return res;
    }

    double norm(const std::vector<double>& v) {
        return std::sqrt(dot(v, v));
    }

    double norm_inf(const std::vector<double>& v) {
        double max_val = 0.0;
        for (double val : v) {
            double abs_val = std::abs(val);
            if (abs_val > max_val) max_val = abs_val;
        }
        return max_val;
    }
}