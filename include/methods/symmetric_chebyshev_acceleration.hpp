#pragma once
#include "methods/iterative_result.hpp"
#include "matrix/plot_matrix.hpp"
#include "vector/operations_vector.hpp"
#include <vector>
#include <functional>
#include <cmath>
#include <chrono>

namespace Methods {

double power_method_max_eigenvalue(
    const PlotMatrix& A,
    size_t max_iterations = 100
);

double estimate_spectral_radius(
    const PlotMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x0,
    size_t num_iter = 20
);

IterativeResult chebyshev_acceleration(
    const PlotMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x,
    double eps,
    double lambda_min,
    double lambda_max,
    size_t chebyshev_degree,
    size_t max_iter
);


template<typename StepMethod>
IterativeResult chebyshev_acceleration_symmetric(
    const PlotMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x0,
    double eps,
    double rho,
    size_t max_iter,
    StepMethod step_method
);


inline std::vector<double> jacobi_step(
    const PlotMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x
) {
    size_t n = A.rows();
    std::vector<double> x_new(n);
    
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += A.get(i, j) * x[j];
            }
        }
        x_new[i] = (b[i] - sum) / A.get(i, i);
    }
    
    return x_new;
}

inline std::vector<double> gauss_seidel_step(
    const PlotMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x
) {
    size_t n = A.rows();
    std::vector<double> x_new = x;
    
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += A.get(i, j) * x_new[j];
            }
        }
        x_new[i] = (b[i] - sum) / A.get(i, i);
    }
    
    return x_new;
}

inline std::vector<double> symmetric_gauss_seidel_step(
    const PlotMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x
) {
    size_t n = A.rows();
    std::vector<double> x_half(n);
    std::vector<double> x_new = x;
    
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += A.get(i, j) * x_new[j];
            }
        }
        x_half[i] = (b[i] - sum) / A.get(i, i);
    }
    
    for (size_t i = n; i-- > 0; ) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += A.get(i, j) * x_half[j];
            }
        }
        x_new[i] = (b[i] - sum) / A.get(i, i);
    }
    
    return x_new;
}

template<typename StepMethod>
IterativeResult chebyshev_acceleration_symmetric(
    const PlotMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x0,
    double eps,
    double rho,
    size_t max_iter,
    StepMethod step_method
) {
    IterativeResult res{ x0, 0, 0.0, 0.0, false, "Chebyshev Symmetric" };
    auto t0 = std::chrono::high_resolution_clock::now();
    
    size_t n = A.rows();
    std::vector<double> y_prev = x0;
    std::vector<double> y_curr = step_method(A, b, x0);
    std::vector<double> y_next(n);
    
    double omega = 1.0;
    double omega_next = 2.0 / (2.0 - rho * rho);
    
    res.iterations = 1;
    
    for (size_t k = 1; k < max_iter; ++k) {
        std::vector<double> step_result = step_method(A, b, y_curr);
        std::vector<double> diff = Vec::sub(step_result, y_prev);
        Vec::mul_num(diff, omega_next);
        y_next = Vec::add(diff, y_prev);
        
        std::vector<double> Ax = A.mul_vec(y_next);
        double nr = Vec::norm(Vec::sub(Ax, b));
        
        res.iterations = k + 1;
        res.residual = nr;
        
        if (nr < eps) {
            res.converged = true;
            break;
        }
        
        y_prev = std::move(y_curr);
        y_curr = std::move(y_next);
        
        omega = omega_next;
        omega_next = 1.0 / (1.0 - (rho * rho * omega) / 4.0);
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    res.x = y_curr;
    
    return res;
}

}  // namespace Methods