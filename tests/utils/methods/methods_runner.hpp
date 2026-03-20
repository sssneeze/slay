#pragma once

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <algorithm>
#include <fstream>

#include "matrix/plot_matrix.hpp"
#include "methods/simple_iteration.hpp"
#include "methods/jacobi.hpp"
#include "methods/gauss_seidel.hpp"
#include "methods/chebyshev_acceleration.hpp"
#include "vector/operations_vector.hpp"

namespace slau {

struct MethodResult {
    std::string name;
    size_t iterations;
    double time_ms;
    double residual;
    bool converged;
};

class MethodsRunner {
public:
    static void run() {
        const size_t N = 200;
        const double eps = 1e-6;
        const double tau = 0.005;
        const size_t max_iter = 10000;

        std::cout << "N = " << N << ", eps = " << eps << ", tau = " << tau << "\n\n";

        auto A = generate_matrix(N);
        auto x_true = generate_vector(N);
        auto b = A.mul_vec(x_true);
        auto x0 = std::vector<double>(N, 0.0);

        std::vector<MethodResult> results;




        auto r1 = Methods::simple_iteration(A, b, x0, tau, eps, max_iter);
        results.push_back({r1.method_name, r1.iterations, r1.time_ms, r1.residual, r1.converged});
        std::cout << "MPI: " << r1.iterations << " iters, " << r1.time_ms << " ms\n";

  
        

        auto r2 = Methods::jacobi(A, b, x0, eps, max_iter);
        results.push_back({r2.method_name, r2.iterations, r2.time_ms, r2.residual, r2.converged});
        std::cout << "Jacobi: " << r2.iterations << " iters, " << r2.time_ms << " ms\n";


        

        auto r3 = Methods::gauss_seidel(A, b, x0, eps, max_iter);
        results.push_back({r3.method_name, r3.iterations, r3.time_ms, r3.residual, r3.converged});
        std::cout << "Gauss-Seidel: " << r3.iterations << " iters, " << r3.time_ms << " ms\n";

        


        double lambda_max = Methods::power_method_max_eigenvalue(A, 100);
        double lambda_min = lambda_max / 10.0;
        auto r4 = Methods::chebyshev_acceleration(A, b, x0, eps, lambda_min, lambda_max, 10, max_iter);
        results.push_back({r4.method_name, r4.iterations, r4.time_ms, r4.residual, r4.converged});
        std::cout << "Chebyshev: " << r4.iterations << " iters, " << r4.time_ms << " ms\n";


        save_results(results, "methods_results.csv");

        auto best = *std::min_element(results.begin(), results.end(),
            [](const MethodResult& a, const MethodResult& b) {
                return a.iterations < b.iterations;
            });
        std::cout << "\nBest: " << best.name << " (" << best.iterations << " iterations)\n";
    }

private:
    static PlotMatrix generate_matrix(size_t n) {
        PlotMatrix A(n, n);
        std::mt19937 gen(42);
        std::uniform_real_distribution<> dis(0.0, 0.5);
        
        for (size_t i = 0; i < n; ++i) {
            double diag = 0.0;
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    A.set(i, j, dis(gen));
                    diag += A.get(i, j);
                }
            }
            A.set(i, i, diag + 1.0);
        }
        return A;
    }

    static std::vector<double> generate_vector(size_t n) {
        std::vector<double> v(n);
        std::mt19937 gen(42);
        std::uniform_real_distribution<> dis(0.0, 10.0);
        for (auto& val : v) val = dis(gen);
        return v;
    }

    static void save_results(const std::vector<MethodResult>& results, const std::string& filename) {
        std::ofstream file(filename);
        file << "method,iterations,time_ms,residual,converged\n";
        for (const auto& r : results) {
            file << r.name << "," << r.iterations << "," << r.time_ms 
                 << "," << r.residual << "," << (r.converged ? "true" : "false") << "\n";
        }
    }
};

} 