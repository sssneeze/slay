#pragma once
#include "matrix/plot_matrix.hpp"
#include <vector>
#include <utility>

namespace Methods {
    double vector_norm(const std::vector<double>& v);
    void householder_transform(const std::vector<double>& v, std::vector<double>& x);
    
    std::pair<PlotMatrix, PlotMatrix> qr_decompose(const PlotMatrix& A);
    
    std::vector<double> solve_with_qr(const PlotMatrix& A, const std::vector<double>& b);
}