#include "methods/simple_iteration.hpp"
#include "vector/operations_vector.hpp"
#include <chrono>

namespace Methods {
    IterativeResult simple_iteration(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double tau,
        double eps,
        size_t max_iter
    ) {
        IterativeResult res{ x, 0, 0.0, 0.0, false, "Simple Iteration" };
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();

        for (size_t k = 0; k < max_iter; ++k) {
            std::vector<double> Ax = A.mul_vec(x);
            std::vector<double> r = Vec::sub(Ax, b);
            double nr = Vec::norm(r);
            
            if (nr < eps) {
                res.converged = true;
                res.iterations = k;
                res.residual = nr;
                break;
            }
            
            for (size_t i = 0; i < n; ++i) {
                x[i] -= tau * r[i];
            }
            
            res.iterations = k + 1;
            res.residual = nr;
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        res.x = x;
        return res;
    }
}