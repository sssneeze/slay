#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Methods {
    std::vector<double> thomas_algorithm(const std::vector<double>& main_d,
                                         const std::vector<double>& lower_d,
                                         const std::vector<double>& upper_d,
                                         const std::vector<double>& rhs,
                                         int N);
}