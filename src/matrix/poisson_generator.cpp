#include "matrix/poisson_generator.hpp"
#include <cmath>

namespace Poisson {
    PlotMatrix generate(size_t nx, size_t ny, double h) {
        size_t N = nx * ny;
        PlotMatrix A(N, N);
        
        double diag = 4.0 / (h * h);
        double off_diag = -1.0 / (h * h);
        
        for (size_t i = 0; i < ny; ++i) {
            for (size_t j = 0; j < nx; ++j) {
                size_t idx = index(i, j, nx);
                
                A.set(idx, idx, diag);
                
                if (j > 0) {
                    A.set(idx, index(i, j - 1, nx), off_diag);
                }
                if (j < nx - 1) {
                    A.set(idx, index(i, j + 1, nx), off_diag);
                }
                if (i > 0) {
                    A.set(idx, index(i - 1, j, nx), off_diag);
                }
                if (i < ny - 1) {
                    A.set(idx, index(i + 1, j, nx), off_diag);
                }
            }
        }
        
        return A;
    }
    
    std::vector<double> generate_rhs(size_t nx, size_t ny, double h) {
        size_t N = nx * ny;
        std::vector<double> b(N, 1.0);
        return b;
    }
}