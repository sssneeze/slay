#include "methods/householder_alg.hpp"
#include <cmath>

namespace Methods {
    double vector_norm(const std::vector<double>& v) {
        double sum = 0.0;
        for (double val : v) {
            sum += val * val;
        }
        return std::sqrt(sum);
    }

    void householder_transform(const std::vector<double>& v, std::vector<double>& x) {
        double mul_vx = 0.0;
        double mul_vv = 0.0;
        
        for (size_t i = 0; i < v.size(); i++) {
            mul_vx += v[i] * x[i];
            mul_vv += v[i] * v[i];
        }
        
        double coef = 2.0 * mul_vx / mul_vv;
        
        for (size_t i = 0; i < v.size(); i++) {
            x[i] = x[i] - coef * v[i];
        }
    }

    std::pair<PlotMatrix, PlotMatrix> qr_decompose(const PlotMatrix& A) {
        size_t m = A.rows();
        size_t n = A.cols();

        PlotMatrix R(m, n);
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                R.set(i, j, A.get(i, j));
            }
        }

        PlotMatrix Q(m, m);
        for (size_t i = 0; i < m; ++i) {
            Q.set(i, i, 1.0);
        }

        for (size_t k = 0; k < m - 1 && k < n; k++) {
            size_t sub_col_size = m - k;
            std::vector<double> x(sub_col_size);
            
            for (size_t i = 0; i < sub_col_size; i++) {
                x[i] = R.get(k + i, k);
            }

            double x_mod = vector_norm(x);
            std::vector<double> v = x;
            v[0] += x_mod;

            // Применяем к столбцам R
            for (size_t i = k; i < n; i++) {
                std::vector<double> col_v(sub_col_size);
                for (size_t j = 0; j < sub_col_size; j++) {
                    col_v[j] = R.get(k + j, i);
                }

                householder_transform(v, col_v);

                for (size_t j = 0; j < sub_col_size; j++) {
                    R.set(k + j, i, col_v[j]);
                }
            }
            
            // Применяем к строкам Q
            for (size_t i = 0; i < m; i++) {
                std::vector<double> row_v(sub_col_size);
                for (size_t j = 0; j < sub_col_size; j++) {
                    row_v[j] = Q.get(i, j + k);
                }

                householder_transform(v, row_v);

                for (size_t j = 0; j < sub_col_size; ++j) {
                    Q.set(i, k + j, row_v[j]);
                }
            }
        }
        
        return {Q, R};
    }

    std::vector<double> solve_with_qr(const PlotMatrix& A, const std::vector<double>& b) {
        size_t n = A.cols();
        auto qr_result = qr_decompose(A);
        PlotMatrix Q = qr_result.first;
        PlotMatrix R = qr_result.second;

        // Вычисляем Q^T * b
        std::vector<double> Qt_b(n, 0.0);
        for (size_t i = 0; i < n; i++) {
            double sum = 0.0;
            for (size_t j = 0; j < n; j++) {
                sum += Q.get(j, i) * b[j];
            }
            Qt_b[i] = sum;
        }

        // Обратная подстановка для Rx = Q^T * b
        std::vector<double> x(n, 0.0);
        for (int i = static_cast<int>(n) - 1; i >= 0; i--) {
            double val = Qt_b[i];
            for (size_t j = i + 1; j < n; j++) {
                val -= R.get(i, j) * x[j];
            }
            x[i] = val / R.get(i, i);
        }

        return x;
    }
}