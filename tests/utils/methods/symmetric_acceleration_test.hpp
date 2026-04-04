// tests/symmetric_acceleration_test.hpp
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

// Структура для хранения данных об итерациях
struct IterationData {
    size_t iteration;
    std::string method;
    double residual;
    double time_ms;
};

// Функция для сохранения данных в CSV
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
    
    // Запуск метода Якоби с сохранением данных по итерациям
    {
        auto res = Methods::jacobi(A, b, x0, EPS, MAX_ITER);
        res.method_name = "Jacobi";
        results.push_back(res);
        
        // Для простоты сохраняем конечную точку
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
    
    // Сохранение данных в CSV
    save_convergence_data(convergence_data, "symmetric_acceleration_results.csv");
    std::cout << "Данные сохранены в symmetric_acceleration_results.csv" << std::endl << std::endl;
    
    // Вывод результатов в консоль
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