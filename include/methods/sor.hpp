#pragma once
#include "matrix/plot_matrix.hpp"
#include "methods/iterative_result.hpp"
#include <vector>

namespace Methods {
    IterativeResult sor(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter,
        double omega = 1.25
    );
}