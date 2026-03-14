#pragma once
#include <vector>
#include <iostream>

class PlotMatrix {
public:
    PlotMatrix(size_t m, size_t n);
    
    double get(size_t i, size_t j) const;
    void set(size_t i, size_t j, double value);
    void print() const;
    
    std::vector<double> mul_vec(const std::vector<double>& v) const;
    
    size_t rows() const;
    size_t cols() const;

private:
    size_t m_;
    size_t n_;
    std::vector<double> data;
};