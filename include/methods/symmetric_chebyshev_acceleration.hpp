#pragma once
#include "iterative_result.hpp"
#include "plot_matrix.hpp"
#include <vector>
#include <functional>
#include <cmath>

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


std::vector<double> jacobi_step(
    const PlotMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x
);


std::vector<double> gauss_seidel_step(
    const PlotMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x
);


std::vector<double> symmetric_gauss_seidel_step(
    const PlotMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x
);

}  // namespace Methods