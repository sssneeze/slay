#include "matrix/plot_matrix.hpp"
#include <iomanip>

PlotMatrix::PlotMatrix(size_t m, size_t n) : m_(m), n_(n) {
    data.resize(m * n, 0.0);
}

double PlotMatrix::get(size_t i, size_t j) const {
    return data[i * n_ + j];
}

void PlotMatrix::set(size_t i, size_t j, double value) {
    data[i * n_ + j] = value;
}

void PlotMatrix::print() const {
    for (size_t i = 0; i < m_; i++) {
        for (size_t j = 0; j < n_; j++) {
            std::cout << std::fixed << std::setprecision(3) 
                      << data[i * n_ + j] << ' ';
        }
        std::cout << std::endl;
    }        
}

std::vector<double> PlotMatrix::mul_vec(const std::vector<double>& v) const {
    std::vector<double> res(m_, 0.0);
    for (size_t i = 0; i < m_; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < n_; j++) {
            sum += data[i * n_ + j] * v[j];
        }
        res[i] = sum;
    }
    return res;
}

size_t PlotMatrix::rows() const { return m_; }
size_t PlotMatrix::cols() const { return n_; }