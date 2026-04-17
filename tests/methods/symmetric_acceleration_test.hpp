#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "methods/iterative_result.hpp"
#include "matrix/plot_matrix.hpp"
#include "methods/jacobi.hpp"
#include "methods/gauss_seidel.hpp"
#include "methods/symmetric_gauss_seidel.hpp"
#include "methods/symmetric_chebyshev_acceleration.hpp"
#include "methods/sor.hpp"
#include "methods/steepest_descent.hpp"
#include "methods/conjugate_gradient.hpp"
#include "matrix/poisson_generator.hpp"

inline PlotMatrix create_poisson_matrix(size_t n) {
    PlotMatrix A(n, n);
    for (size_t i = 0; i < n; ++i) {
        A.set(i, i, 4.0);
        if (i > 0) A.set(i, i - 1, -1.0);
        if (i < n - 1) A.set(i, i + 1, -1.0);
    }
    return A;
}

inline std::vector<double> create_rhs(size_t n) {
    return std::vector<double>(n, 1.0);
}

inline std::vector<double> create_initial_guess(size_t n) {
    return std::vector<double>(n, 0.0);
}


struct IterationData {
    size_t iteration;
    std::string method;
    double residual;
    double time_ms;
};


inline void save_convergence_data(const std::vector<IterationData>& data, const std::string& filename) {
    std::ofstream file(filename);
    file << "iteration,method,residual,time_ms\n";
    for (const auto& d : data) {
        file << d.iteration << "," << d.method << "," << d.residual << "," << d.time_ms << "\n";
    }
}

inline int run_symmetric_acceleration_test() {
    const size_t N = 100;
    const double EPS = 1e-6;
    const size_t MAX_ITER = 10000;
    
    std::cout << "Сравнение итерационных методов" << std::endl;
    std::cout << "Размер матрицы: " << N << "x" << N << std::endl;
    std::cout << "Точность: " << EPS << std::endl << std::endl;
    
    PlotMatrix A = create_poisson_matrix(N);
    std::vector<double> b = create_rhs(N);
    std::vector<double> x0 = create_initial_guess(N);
    
    double rho = Methods::estimate_spectral_radius(A, b, x0, 20);
    std::cout << "Оценка спектрального радиуса: " << rho << std::endl << std::endl;
    
    std::vector<Methods::IterativeResult> results;
    std::vector<IterationData> convergence_data;
    
    
    {
        auto res = Methods::jacobi(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Jacobi";
        results.push_back(res);
        
        
        convergence_data.push_back({res.iterations, "Jacobi", res.residual, res.time_ms});
    }
    
    // Запуск метода Гаусса-Зейделя
    {
        auto res = Methods::gauss_seidel(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Gauss-Seidel";
        results.push_back(res);
        
        convergence_data.push_back({res.iterations, "Gauss-Seidel", res.residual, res.time_ms});
    }
    
    // Запуск симметризованного Гаусса-Зейделя
    {
        auto res = Methods::symmetric_gauss_seidel(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Symmetric GS";
        results.push_back(res);
        
        convergence_data.push_back({res.iterations, "Symmetric GS", res.residual, res.time_ms});
    }
    
    // Ускорение Чебышёва для Якоби
    {
        auto res = Methods::chebyshev_acceleration_symmetric(
            A, b, x0, EPS, rho, MAX_ITER, Methods::jacobi_step
        );
        res.method_name = "Chebyshev + Jacobi";
        results.push_back(res);
        
        convergence_data.push_back({res.iterations, "Chebyshev + Jacobi", res.residual, res.time_ms});
    }
    
    // Ускорение Чебышёва для Гаусса-Зейделя
    {
        auto res = Methods::chebyshev_acceleration_symmetric(
            A, b, x0, EPS, rho, MAX_ITER, Methods::gauss_seidel_step
        );
        res.method_name = "Chebyshev + GS";
        results.push_back(res);
        
        convergence_data.push_back({res.iterations, "Chebyshev + GS", res.residual, res.time_ms});
    }
    
    // Ускорение Чебышёва для симметризованного Г-З
    {
        auto res = Methods::chebyshev_acceleration_symmetric(
            A, b, x0, EPS, rho, MAX_ITER, Methods::symmetric_gauss_seidel_step
        );
        res.method_name = "Chebyshev + SGS";
        results.push_back(res);
        
        convergence_data.push_back({res.iterations, "Chebyshev + SGS", res.residual, res.time_ms});
    }
    
    
    save_convergence_data(convergence_data, "symmetric_acceleration_results.csv");
    std::cout << "Данные сохранены в symmetric_acceleration_results.csv" << std::endl << std::endl;
    
    
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
                  << (res.converged ? ":)" : ":(") << std::endl;
    }
    
    std::cout << std::endl << "Ускорение относительно Якоби" << std::endl;
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
//======================================================================================

    inline int run_seminar7_test() {
    const size_t NX = 20, NY = 20;
    const double EPS = 1e-6;
    const size_t MAX_ITER = 5000;
    const double H = 1.0;
    
    std::cout << "Сравнение методов на матрице Пуассона (" << NX << "×" << NY << ")" << std::endl;
    std::cout << "Точность: " << EPS << ", макс. итераций: " << MAX_ITER << std::endl << std::endl;
    
    PlotMatrix A = Poisson::generate(NX, NY, H);
    std::vector<double> b = Poisson::generate_rhs(NX, NY, H);
    std::vector<double> x0(NX * NY, 0.0);
    
    std::vector<Methods::IterativeResult> results;
    
    {
        auto res = Methods::gauss_seidel(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Gauss-Seidel";
        results.push_back(res);
    }
    
    {
        auto res = Methods::sor(A, b, x0, EPS, MAX_ITER, 1.5);
        res.method_name = "SOR (omega=1.5)";
        results.push_back(res);
    }
    
    {
        double rho = Methods::estimate_spectral_radius(A, b, x0, 20);
        auto res = Methods::chebyshev_acceleration_symmetric(
            A, b, x0, EPS, rho, MAX_ITER, Methods::gauss_seidel_step
        );
        res.method_name = "Chebyshev + GS";
        results.push_back(res);
    }
    
    {
        auto res = Methods::steepest_descent(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Steepest Descent";
        results.push_back(res);
    }
    
    {
        auto res = Methods::conjugate_gradient(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Conjugate Gradient";
        results.push_back(res);
    }
    
    std::cout << std::left << std::setw(25) << "Метод" 
              << std::setw(12) << "Итерации" 
              << std::setw(15) << "Время (мс)" 
              << std::setw(18) << "Невязка" 
              << "Сходимость" << std::endl;
    std::cout << std::string(90, '-') << std::endl;
    
    for (const auto& res : results) {
        std::cout << std::left << std::setw(25) << res.method_name
                  << std::setw(12) << res.iterations
                  << std::setw(15) << std::fixed << std::setprecision(2) << res.time_ms
                  << std::setw(18) << std::scientific << std::setprecision(2) << res.residual
                  << (res.converged ? ":)" : ":(") << std::endl;
    }
    
    std::cout << "\nУскорение относительно Гаусс-Зейделя:" << std::endl;
    double base_iter = results[0].iterations;
    double base_time = results[0].time_ms;
    
    for (size_t i = 1; i < results.size(); ++i) {
        if (results[i].iterations > 0 && results[i].time_ms > 0) {
            double speedup_iter = static_cast<double>(base_iter) / results[i].iterations;
            double speedup_time = base_time / results[i].time_ms;
            std::cout << "  " << results[i].method_name << ": "
                      << std::fixed << std::setprecision(2)
                      << speedup_iter << "x итер., "
                      << speedup_time << "x время" << std::endl;
        }
    }
    
    return 0;

}