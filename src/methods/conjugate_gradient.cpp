#include "methods/conjugate_gradient.hpp"
#include "vector/operations_vector.hpp"
#include <chrono>

namespace Methods {
    IterativeResult conjugate_gradient(
        const PlotMatrix& A,
        const std::vector<double>& b,
        std::vector<double> x,
        double eps,
        size_t max_iter
    ) {
        IterativeResult res{ x, 0, 0.0, 0.0, false, "Conjugate Gradient" };
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();

        std::vector<double> Ax = A.mul_vec(x);
        std::vector<double> r = Vec::sub(b, Ax);
        std::vector<double> d = r;
        
        double r_dot_r = Vec::dot(r, r);

        for (size_t k = 0; k < max_iter; ++k) {
            std::vector<double> Ad = A.mul_vec(d);
            double d_dot_Ad = Vec::dot(d, Ad);
            
            if (std::abs(d_dot_Ad) < 1e-15) break;
            
            double alpha = r_dot_r / d_dot_Ad;
            
            for (size_t i = 0; i < n; ++i) {
                x[i] += alpha * d[i];
            }
            
            std::vector<double> r_new(n);
            for (size_t i = 0; i < n; ++i) {
                r_new[i] = r[i] - alpha * Ad[i];
            }
            
            double nr = Vec::norm(r_new);
            res.iterations = k + 1;
            res.residual = nr;
            
            if (nr < eps) {
                res.converged = true;
                break;
            }
            
            double r_new_dot_r_new = Vec::dot(r_new, r_new);
            double beta = r_new_dot_r_new / r_dot_r;
            
            for (size_t i = 0; i < n; ++i) {
                d[i] = r_new[i] + beta * d[i];
            }
            
            r = std::move(r_new);
            r_dot_r = r_new_dot_r_new;
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        res.x = x;
        return res;
    }
}