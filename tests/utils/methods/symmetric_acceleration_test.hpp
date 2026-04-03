#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "iterative_result.hpp"
#include "matrix/plot_matrix.hpp"
#include "methods/jacobi.hpp"
#include "methods/gauss_seidel.hpp"
#include "methods/symmetric_gauss_seidel.hpp"
#include "methods/chebyshev_acceleration.hpp"

PlotMatrix create_poisson_matrix(size_t n) {
    PlotMatrix A(n, n);
    
    for (size_t i = 0; i < n; ++i) {
        A.set(i, i, 4.0);
        if (i > 0) A.set(i, i - 1, -1.0);
        if (i < n - 1) A.set(i, i + 1, -1.0);
    }
    
    return A;
}

std::vector<double> create_rhs(size_t n) {
    return std::vector<double>(n, 1.0);
}

std::vector<double> create_initial_guess(size_t n) {
    return std::vector<double>(n, 0.0);
}

int run_main() {
    const size_t N = 100;
    const double EPS = 1e-6;
    const size_t MAX_ITER = 10000;
    
    std::cout << "=== Сравнение итерационных методов ===" << std::endl;
    std::cout << "Размер матрицы: " << N << "x" << N << std::endl;
    std::cout << "Точность: " << EPS << std::endl;
    std::cout << std::endl;
    
    PlotMatrix A = create_poisson_matrix(N);
    std::vector<double> b = create_rhs(N);
    std::vector<double> x0 = create_initial_guess(N);
    
    double rho = Methods::estimate_spectral_radius(A, b, x0, 20);
    std::cout << "Оценка спектрального радиуса: " << rho << std::endl;
    std::cout << std::endl;
    
    std::vector<Methods::IterativeResult> results;
    
    results.push_back(Methods::jacobi(A, b, x0, EPS, MAX_ITER));
    results.push_back(Methods::gauss_seidel(A, b, x0, EPS, MAX_ITER));
    results.push_back(Methods::symmetric_gauss_seidel(A, b, x0, EPS, MAX_ITER));
    
    auto cheb_jacobi = Methods::chebyshev_acceleration_symmetric(
        A, b, x0, EPS, rho, MAX_ITER, Methods::jacobi_step
    );
    cheb_jacobi.method_name = "Chebyshev + Jacobi";
    results.push_back(cheb_jacobi);
    
    auto cheb_gs = Methods::chebyshev_acceleration_symmetric(
        A, b, x0, EPS, rho, MAX_ITER, Methods::gauss_seidel_step
    );
    cheb_gs.method_name = "Chebyshev + GS";
    results.push_back(cheb_gs);
    
    auto cheb_sgs = Methods::chebyshev_acceleration_symmetric(
        A, b, x0, EPS, rho, MAX_ITER, Methods::symmetric_gauss_seidel_step
    );
    cheb_sgs.method_name = "Chebyshev + SGS";
    results.push_back(cheb_sgs);
    
    std::cout << std::left << std::setw(25) << "Метод" 
              << std::setw(15) << "Итерации" 
              << std::setw(15) << "Время (мс)" 
              << std::setw(15) << "Невязка" 
              << "Сходимость" << std::endl;
    std::cout << std::string(85, '-') << std::endl;
    
    for (const auto& res : results) {
        std::cout << std::left << std::setw(25) << res.method_name
                  << std::setw(15) << res.iterations
                  << std::setw(15) << std::fixed << std::setprecision(2) << res.time_ms
                  << std::setw(15) << std::scientific << std::setprecision(2) << res.residual
                  << (res.converged ? "✓" : "✗") << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "=== Ускорение относительно Якоби ===" << std::endl;
    double base_iter = results[0].iterations;
    double base_time = results[0].time_ms;
    
    for (size_t i = 1; i < results.size(); ++i) {
        double speedup_iter = base_iter / results[i].iterations;
        double speedup_time = base_time / results[i].time_ms;
        std::cout << results[i].method_name << ": "
                  << std::fixed << std::setprecision(2)
                  << speedup_iter << "x (итерации), "
                  << speedup_time << "x (время)" << std::endl;
    }
    
    return 0;
}