#pragma once
#include "matrix/plot_matrix.hpp"
#include "methods/iterative_result.hpp"
#include <vector>

namespace Methods {
    IterativeResult steepest_descent(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter
    );
}