// include/matrix/csr_matrix.hpp
#pragma once
#include <vector>
#include <tuple>
#include <iostream>

class CSRMatrix {
public:
    CSRMatrix(size_t m, size_t n, 
              const std::vector<std::tuple<int, int, double>>& dok);

    std::vector<double> mul_vec(const std::vector<double>& v) const;
    void print() const;
    

    size_t rows() const { return m_; }
    size_t cols() const { return n_; }
    

    const std::vector<double>& get_values() const { return values_; }
    const std::vector<int>& get_cols() const { return cols_; }
    const std::vector<int>& get_rows() const { return rows_; }

private:
    size_t m_;  // число строк
    size_t n_;  // число столбцов
    std::vector<double> values_;  // ненулевые значения
    std::vector<int> cols_;       // номера столбцов
    std::vector<int> rows_;       // указатели на строки
};