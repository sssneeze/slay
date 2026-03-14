#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <cmath>

#include "matrix/plot_matrix.hpp"
#include "methods/simple_iteration.hpp"
#include "methods/jacobi.hpp"
#include "methods/gauss_seidel.hpp"
#include "vector/operations_vector.hpp"

namespace slau {

struct MethodStats {
    std::string name;
    size_t iterations;
    double time_ms;
    double residual;
    bool converged;
    std::vector<double> solution;
};

class MethodsRunner {
public:
    static void run() {

        const size_t N = 200;
        const double eps = 1e-6;
        const double tau = 0.005;
        const size_t max_iter = 10000;

        std::cout << "Параметры теста:" << std::endl;
        std::cout << "Размер матрицы: " << N << "×" << N << std::endl;
        std::cout << "Точность: " << eps << std::endl;
        std::cout << "Параметр tau: " << tau << std::endl;
        std::cout << "Макс. итераций: " << max_iter << std::endl;
        std::cout << std::endl;


        auto A = generate_diagonally_dominant_matrix(N);
        auto x_true = generate_random_vector(N, 0.0, 10.0);
        auto b = A.mul_vec(x_true);
        auto x0 = std::vector<double>(N, 0.0);

        std::vector<MethodStats> results;

        //метод мпи
        results.push_back(run_simple_iteration(A, b, x0, tau, eps, max_iter));
        std::cout << results.back().iterations << " итераций" << std::endl;

        // метод якоби
        results.push_back(run_jacobi(A, b, x0, eps, max_iter));
        std::cout << results.back().iterations << " итераций" << std::endl;

        // метод гаусса-зейделя
        results.push_back(run_gauss_seidel(A, b, x0, eps, max_iter));
        std::cout << results.back().iterations << " итераций" << std::endl;
        std::cout << std::endl;

        print_results_table(results);
    }

private:
    static PlotMatrix generate_diagonally_dominant_matrix(size_t n) {
        PlotMatrix A(n, n);
        std::mt19937 gen(42); 
        std::uniform_real_distribution<> dis(0.0, 0.5);
        
        for (size_t i = 0; i < n; ++i) {
            double diag = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    double val = dis(gen);
                    A.set(i, j, val);
                    diag += val;
                }
            }
            A.set(i, i, diag + 1.0);
        }
        return A;
    }

    static std::vector<double> generate_random_vector(size_t n, double min_val, double max_val) {
        std::vector<double> v(n);
        std::mt19937 gen(42);
        std::uniform_real_distribution<> dis(min_val, max_val);
        for (auto& val : v) val = dis(gen);
        return v;
    }

    static MethodStats run_simple_iteration(const PlotMatrix& A, const std::vector<double>& b,
                                           std::vector<double> x, double tau, double eps, 
                                           size_t max_iter) {
        MethodStats stats{"Simple Iteration", 0, 0.0, 0.0, false, {}};
        
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();

        for (size_t k = 0; k < max_iter; ++k) {
            auto Ax = A.mul_vec(x);
            auto r = Vec::sub(Ax, b);
            double nr = Vec::norm(r);
            
            if (nr < eps) {
                stats.converged = true;
                stats.iterations = k;
                stats.residual = nr;
                break;
            }
            
            for (size_t i = 0; i < n; ++i) {
                x[i] -= tau * r[i];
            }
            
            stats.iterations = k + 1;
            stats.residual = nr;
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        stats.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        stats.solution = x;
        return stats;
    }

    static MethodStats run_jacobi(const PlotMatrix& A, const std::vector<double>& b,
                                 std::vector<double> x, double eps, size_t max_iter) {
        MethodStats stats{"Jacobi", 0, 0.0, 0.0, false, {}};
        
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();
        std::vector<double> x_new(n);

        for (size_t k = 0; k < max_iter; ++k) {
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    if (i != j) {
                        sum += A.get(i, j) * x[j];
                    }
                }
                x_new[i] = (b[i] - sum) / A.get(i, i);
            }
            
            auto Ax = A.mul_vec(x_new);
            double nr = Vec::norm(Vec::sub(Ax, b));
            x = x_new;
            stats.iterations = k + 1;
            stats.residual = nr;
            
            if (nr < eps) {
                stats.converged = true;
                break;
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        stats.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        stats.solution = x;
        return stats;
    }

    static MethodStats run_gauss_seidel(const PlotMatrix& A, const std::vector<double>& b,
                                       std::vector<double> x, double eps, size_t max_iter) {
        MethodStats stats{"Gauss-Seidel", 0, 0.0, 0.0, false, {}};
        
        auto t0 = std::chrono::high_resolution_clock::now();
        size_t n = A.rows();

        for (size_t k = 0; k < max_iter; ++k) {
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    if (i != j) {
                        sum += A.get(i, j) * x[j];
                    }
                }
                x[i] = (b[i] - sum) / A.get(i, i);
            }
            
            auto Ax = A.mul_vec(x);
            double nr = Vec::norm(Vec::sub(Ax, b));
            stats.iterations = k + 1;
            stats.residual = nr;
            
            if (nr < eps) {
                stats.converged = true;
                break;
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        stats.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        stats.solution = x;
        return stats;
    }

    static void print_results_table(const std::vector<MethodStats>& results) {
        for (const auto& stat : results) {
            std::cout << std::left << std::setw(24) << stat.name << "  " 
                      << std::right << std::setw(10) << stat.iterations << "  "
                      << std::setw(11) << std::fixed << std::setprecision(2) << stat.time_ms << "  "
                      << std::scientific << std::setprecision(2) << std::setw(13) << stat.residual << " \n";
        }
        std::cout << std::endl;
    }
};

} 