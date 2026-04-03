#include "symmetric_chebyshev_acceleration.hpp"
#include "operations_vector.hpp"
#include <chrono>
#include <random>
#include <utility>
#include <numbers>

namespace Methods {

double power_method_max_eigenvalue(
    const PlotMatrix& A,
    size_t max_iterations
) {
    size_t n = A.rows();
    
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::vector<double> r(n);
    for (auto& val : r) val = dis(gen);
    
    double lambda_max = 0.0;
    
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        auto Ar = A.mul_vec(r);
        
        double norm = Vec::norm(Ar);
        if (norm < 1e-15) break;
        
        double numerator = 0.0, denominator = 0.0;
        for (size_t i = 0; i < n; ++i) {
            numerator += r[i] * Ar[i];
            denominator += r[i] * r[i];
        }
        
        lambda_max = numerator / denominator;
        
        for (size_t i = 0; i < n; ++i) {
            r[i] = Ar[i] / norm;
        }
    }
    
    return lambda_max;
}

double estimate_spectral_radius(
    const PlotMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x0,
    size_t num_iter
) {
    size_t n = A.rows();
    std::vector<double> x_prev = x0;
    std::vector<double> x_curr = x0;
    std::vector<double> x_next(n);
    
    for (size_t k = 0; k < num_iter; ++k) {
        x_prev = x_curr;
        for (size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A.get(i, j) * x_curr[j];
                }
            }
            x_next[i] = (b[i] - sum) / A.get(i, i);
        }
        x_curr = x_next;
    }
    
    std::vector<double> x_next2(n);
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += A.get(i, j) * x_curr[j];
            }
        }
        x_next2[i] = (b[i] - sum) / A.get(i, i);
    }
    
    std::vector<double> diff1 = Vec::sub(x_curr, x_prev);
    std::vector<double> diff2 = Vec::sub(x_next2, x_curr);
    
    double norm1 = Vec::norm(diff1);
    double norm2 = Vec::norm(diff2);
    
    if (norm1 < 1e-12) return 0.5;
    
    double rho = norm2 / norm1;
    return std::min(rho, 0.99);
}

IterativeResult chebyshev_acceleration(
    const PlotMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x,
    double eps,
    double lambda_min,
    double lambda_max,
    size_t chebyshev_degree,
    size_t max_iter
) {
    IterativeResult res{ x, 0, 0.0, 0.0, false, "Chebyshev Acceleration" };
    auto t0 = std::chrono::high_resolution_clock::now();
    size_t n = A.rows();
    
    double center = (lambda_min + lambda_max) / 2.0;
    double half_width = (lambda_max - lambda_min) / 2.0;
    
    size_t total_iter = 0;
    bool converged = false;
    
    while (total_iter < max_iter && !converged) {
        size_t steps_in_cycle = std::min(chebyshev_degree, max_iter - total_iter);
        
        std::vector<double> taus(steps_in_cycle);
        std::vector<size_t> permutation(steps_in_cycle);
        
        auto build_permutation = [](size_t n) {
            std::vector<size_t> perm(n);
            if (n == 1) {
                perm[0] = 0;
                return perm;
            }
            
            size_t m = 1;
            while (m < n) m *= 2;
            
            std::vector<size_t> temp(m);
            temp[0] = 0;
            temp[1] = 1;
            
            size_t current_size = 2;
            while (current_size < m) {
                std::vector<size_t> new_temp(current_size * 2);
                for (size_t i = 0; i < current_size; ++i) {
                    new_temp[2*i] = temp[i];
                    new_temp[2*i + 1] = current_size * 2 - 1 - temp[i];
                }
                temp = std::move(new_temp);
                current_size *= 2;
            }
            
            for (size_t i = 0; i < n; ++i) {
                perm[i] = temp[i];
            }
            return perm;
        };
        
        permutation = build_permutation(steps_in_cycle);
        
        for (size_t k = 0; k < steps_in_cycle; ++k) {
            size_t idx = permutation[k];
            double t_k = std::cos(std::numbers::pi * (2.0 * idx + 1.0) / (2.0 * steps_in_cycle));
            double lambda_k = center + half_width * t_k;
            taus[k] = 1.0 / lambda_k;
        }
        
        for (size_t k = 0; k < steps_in_cycle; ++k) {
            double tau = taus[k];
            auto Ax = A.mul_vec(x);
            auto r = Vec::sub(Ax, b);
            
            for (size_t i = 0; i < n; ++i) {
                x[i] -= tau * r[i];
            }
            
            total_iter++;
            double nr = Vec::norm(r);
            if (nr < eps) {
                converged = true;
                res.residual = nr;
                break;
            }
        }
        
        if (converged) break;
        
        auto Ax = A.mul_vec(x);
        double nr = Vec::norm(Vec::sub(Ax, b));
        if (nr < eps) {
            converged = true;
            res.residual = nr;
        }
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    res.iterations = total_iter;
    res.converged = converged;
    res.x = x;
    
    return res;
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

std::vector<double> jacobi_step(
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

std::vector<double> gauss_seidel_step(
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

std::vector<double> symmetric_gauss_seidel_step(
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
    
    for (int i = n - 1; i >= 0; --i) {
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

}