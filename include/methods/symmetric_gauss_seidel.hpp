#pragma once
#include "methods/iterative_result.hpp"
#include "../matrix/plot_matrix.hpp"
#include <vector>
#include <chrono>

namespace Methods {
    IterativeResult symmetric_gauss_seidel(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter
    );
}