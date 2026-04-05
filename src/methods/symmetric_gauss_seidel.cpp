#include "methods/symmetric_gauss_seidel.hpp"
#include "vector/operations_vector.hpp"

namespace Methods {
    IterativeResult symmetric_gauss_seidel(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter
    ) {
        IterativeResult res{ x, 0, 0.0, 0.0, false, "Symmetric GS" };
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();
        std::vector<double> x_half(n);

        for (size_t k = 0; k < max_iter; ++k) {
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    if (i != j) {
                        sum += A.get(i, j) * x[j];
                    }
                }
                x_half[i] = (b[i] - sum) / A.get(i, i);
            }
            
            for (int i = n - 1; i >= 0; --i) {
                double sum = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    if (i != j) {
                        sum += A.get(i, j) * x_half[j];
                    }
                }
                x[i] = (b[i] - sum) / A.get(i, i);
            }
            
            std::vector<double> Ax = A.mul_vec(x);
            double nr = Vec::norm(Vec::sub(Ax, b));
            res.iterations = k + 1;
            res.residual = nr;
            
            if (nr < eps) {
                res.converged = true;
                break;
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        res.x = x;
        return res;
    }
}