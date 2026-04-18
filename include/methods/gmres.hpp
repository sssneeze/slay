#pragma once
#include "matrix/csr_matrix.hpp"
#include "methods/iterative_result.hpp"
#include <vector>

namespace Methods {
    IterativeResult gmres(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter,
        size_t restart = 30
    );
}