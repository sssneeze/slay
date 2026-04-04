#pragma once
#include <vector>
#include <string>

namespace Methods {
    struct IterativeResult {
        std::vector<double> x;          
        size_t iterations;             
        double time_ms;                 
        double residual;               
        bool converged;               
        std::string method_name;      
    };
}