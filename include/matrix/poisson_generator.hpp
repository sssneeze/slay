#pragma once
#include "matrix/plot_matrix.hpp"
#include <vector>

namespace Poisson {
    PlotMatrix generate(size_t nx, size_t ny, double h = 1.0);
    std::vector<double> generate_rhs(size_t nx, size_t ny, double h = 1.0);
    
    inline size_t index(size_t i, size_t j, size_t nx) {
        return i * nx + j;
    }
}