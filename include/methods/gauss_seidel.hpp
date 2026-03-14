#pragma once
#include "matrix/plot_matrix.hpp"
#include "methods/iterative_result.hpp"
#include <vector>

namespace Methods {
    IterativeResult gauss_seidel(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter = 10000
    );
}