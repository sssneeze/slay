#include "methods/gmres.hpp"
#include "vector/operations_vector.hpp"
#include <chrono>
#include <cmath>
#include <algorithm>

namespace Methods {

using GivensRotation = std::pair<double, double>;

struct HessenbergMatrix {
    size_t max_m;
    std::vector<double> data;
    
    size_t index(size_t i, size_t j) const noexcept {
        if (i > j + 1) return static_cast<size_t>(-1);
        size_t idx = 0;
        for (size_t col = 0; col < j; ++col) {
            idx += std::min(col + 2, max_m);
        }
        return idx + i;
    }
    
    double& operator()(size_t i, size_t j) { return data[index(i, j)]; }
    const double& operator()(size_t i, size_t j) const { return data[index(i, j)]; }
    
    HessenbergMatrix(size_t m) : max_m(m) {
        size_t total = 0;
        for (size_t j = 0; j < m; ++j) {
            total += std::min(j + 2, m);
        }
        data.resize(total, 0.0);
    }
};

inline GivensRotation compute_givens(double a, double b) {
    if (std::abs(b) < 1e-15) return {1.0, 0.0};
    if (std::abs(a) < 1e-15) return {0.0, b > 0 ? 1.0 : -1.0};
    double r = std::hypot(a, b);
    return {a / r, -b / r};
}

inline void apply_givens(const GivensRotation& rot, double& a, double& b) {
    double c = rot.first, s = rot.second;
    double ta = c * a - s * b;
    double tb = s * a + c * b;
    a = ta; b = tb;
}

inline void solve_upper_triangular(
    const HessenbergMatrix& R,
    const std::vector<double>& rhs,
    std::vector<double>& result,
    size_t m
) {
    for (size_t i = m; i-- > 0; ) {
        double sum = rhs[i];
        for (size_t j = i + 1; j < m; ++j) {
            sum -= R(i, j) * result[j];
        }
        double diag = R(i, i);
        result[i] = (std::abs(diag) < 1e-15) ? 0.0 : sum / diag;
    }
}

IterativeResult gmres(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double> x,
    double eps,
    size_t max_iter,
    size_t restart
) {
    IterativeResult res{ x, 0, 0.0, 0.0, false, "GMRES" };
    auto t0 = std::chrono::high_resolution_clock::now();
    
    size_t n = A.rows();
    if (restart == 0) restart = n;
    restart = std::min(restart, n);
    
    std::vector<double> Ax = A.mul_vec(x);
    std::vector<double> r = Vec::sub(b, Ax);
    double beta = Vec::norm(r);
    
    if (beta < eps) {
        res.converged = true;
        res.residual = beta;
        auto t1 = std::chrono::high_resolution_clock::now();
        res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        return res;
    }
    
    size_t total_iter = 0;
    bool converged = false;
    
    while (total_iter < max_iter && !converged) {
        size_t m = std::min(restart, max_iter - total_iter);
        
        std::vector<std::vector<double>> V(m + 1, std::vector<double>(n));
        for (size_t i = 0; i < n; ++i) {
            V[0][i] = r[i] / beta;
        }
        
        HessenbergMatrix H(m);
        std::vector<GivensRotation> rotations;
        std::vector<double> g(m + 1, 0.0);
        g[0] = beta;
        
        for (size_t j = 0; j < m && !converged; ++j) {
            std::vector<double> w = A.mul_vec(V[j]);
            
            for (size_t i = 0; i <= j; ++i) {
                H(i, j) = Vec::dot(V[i], w);
                for (size_t k = 0; k < n; ++k) {
                    w[k] -= H(i, j) * V[i][k];
                }
            }
            
            double h_norm = Vec::norm(w);
            H(j + 1, j) = h_norm;
            
            if (h_norm > 1e-15) {
                for (size_t k = 0; k < n; ++k) {
                    V[j + 1][k] = w[k] / h_norm;
                }
            }
            
            for (size_t i = 0; i < j; ++i) {
                apply_givens(rotations[i], H(i, j), H(i + 1, j));
            }
            
            GivensRotation rot = compute_givens(H(j, j), H(j + 1, j));
            rotations.push_back(rot);
            apply_givens(rot, H(j, j), H(j + 1, j));
            apply_givens(rot, g[j], g[j + 1]);
            
            double residual = std::abs(g[j + 1]);
            total_iter++;
            
            if (residual < eps) {
                std::vector<double> y(j + 1);
                solve_upper_triangular(H, g, y, j + 1);
                
                for (size_t i = 0; i <= j; ++i) {
                    for (size_t k = 0; k < n; ++k) {
                        x[k] += y[i] * V[i][k];
                    }
                }
                
                converged = true;
                res.residual = residual;
                break;
            }
        }
        
        if (!converged) {
            Ax = A.mul_vec(x);
            r = Vec::sub(b, Ax);
            beta = Vec::norm(r);
            
            if (beta < eps) {
                converged = true;
                res.residual = beta;
            }
        }
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    res.iterations = total_iter;
    res.converged = converged;
    res.x = x;
    
    return res;
}

}  // namespace Methods