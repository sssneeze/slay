#pragma once
#include <vector>
#include <string>

struct IterativeResult {
    std::vector<double> x;          // решение
    size_t iterations;              // число итераций
    double time_ms;                 // время в мс
    double residual;                // финальная невязка ||Ax-b||
    bool converged;                 // достигнута ли точность
    std::string method_name;        // имя метода
};