#include "methods/steepest_descent.hpp"
#include "vector/operations_vector.hpp"
#include <chrono>

namespace Methods {
    IterativeResult steepest_descent(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter
    ) {
        IterativeResult res{ x, 0, 0.0, 0.0, false, "Steepest Descent" };
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();

        for (size_t k = 0; k < max_iter; ++k) {
            std::vector<double> Ax = A.mul_vec(x);
            std::vector<double> r = Vec::sub(b, Ax);
            
            std::vector<double> Ar = A.mul_vec(r);
            double r_dot_r = Vec::dot(r, r);
            double r_dot_Ar = Vec::dot(r, Ar);
            
            if (std::abs(r_dot_Ar) < 1e-15) break;
            
            double alpha = r_dot_r / r_dot_Ar;
            
            for (size_t i = 0; i < n; ++i) {
                x[i] += alpha * r[i];
            }
            
            double nr = Vec::norm(r);
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