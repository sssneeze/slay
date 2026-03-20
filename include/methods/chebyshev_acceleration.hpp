#pragma once
#include "matrix/plot_matrix.hpp"
#include "methods/iterative_result.hpp"
#include "vector/operations_vector.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace Methods {
    IterativeResult chebyshev_acceleration(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        double lambda_min,
        double lambda_max,
        size_t chebyshev_degree = 10,
        size_t max_iter = 10000
    );

    double power_method_max_eigenvalue(
        const PlotMatrix& A,
        size_t max_iterations = 100
    );
}