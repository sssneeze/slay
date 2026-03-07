#include <iostream>
#include <vector>
#include <cmath>
#include <utility>


class PlotMatrix {
public:
    PlotMatrix(size_t m, size_t n) : m_(m), n_(n) {
        data.resize(m * n, 0.0);
    }

    double get(size_t i, size_t j) const {
        return data[i * n_ + j];
    }

    void set(size_t i, size_t j, double value) {
        data[i * n_ + j] = value;
    }

    void print() const {
        for (size_t i = 0; i < m_; i++) {
            for (size_t j = 0; j < n_; j++) {
                std::cout << data[i * n_ + j] << ' ';
            }
            std::cout << std::endl;
        }        
    }

    
    std::vector<double> mul_vec(const std::vector<double>& v) const {
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

    size_t rows() const { return m_; }
    size_t cols() const { return n_; }

private:
    size_t m_; 
    size_t n_; 
    std::vector<double> data;
};


double vector_mod(const std::vector<double>& v) {
    double sum = 0.0;
    for (size_t i = 0; i < v.size(); i++) {
        sum += v[i] * v[i];
    }
    return std::sqrt(sum);
}


void householder_alg(const std::vector<double>& v, std::vector<double>& x) {
    
    double mul_vx = 0.0;
    double mul_vv = 0.0;
    for (size_t i = 0; i < v.size(); i++) {
        mul_vx += v[i] * x[i];
        mul_vv += v[i] * v[i];
    }
    double coef = 2.0 * mul_vx/mul_vv; 
    for (size_t i = 0; i < v.size(); i++) {
        x[i] = x[i] - coef * v[i];
    }
}


std::pair<PlotMatrix, PlotMatrix> QR_alg(const PlotMatrix& A) {
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

        double x_mod = vector_mod(x);
        std::vector<double> v = x;
        v[0] += x_mod;

        for (size_t i = k; i < n; i++) {
            std::vector<double> col_v(sub_col_size);
            for (size_t j = 0; j < sub_col_size; j++) {
                col_v[j] = R.get(k + j, i); 
            }

            householder_alg(v, col_v);

            for (size_t j = 0; j < sub_col_size; j++) {
                R.set(k + j, i, col_v[j]); 
            }
        }
        
        for (size_t i = 0; i < m; i++) {
            std::vector<double> row_v(sub_col_size);
            for (size_t j = 0; j < sub_col_size; j++) {
                row_v[j] = Q.get(i, j + k);
            }

            householder_alg(v, row_v);

            for (size_t j = 0; j < sub_col_size; ++j) {
                Q.set(i, k + j, row_v[j]);
            }
        }
    }
    return {Q, R};
}


std::vector<double> solve_slay(const PlotMatrix& A, const std::vector<double>& b) {
    size_t n = A.cols();
    auto qr_result = QR_alg(A);
    PlotMatrix Q = qr_result.first;
    PlotMatrix R = qr_result.second;

    std::vector<double> Qt_b(n, 0.0);
    for (size_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < n; j++) {
            sum += Q.get(j, i) * b[j];
        }
        Qt_b[i] = sum;
    }

    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double val = Qt_b[i];

        for (size_t j = i + 1; j < n; j++) {
            val = val - R.get(i, j) * x[j];
        }

        x[i] = val / R.get(i, i);
    }

    return x;
}


int main() {
    PlotMatrix A(3, 3);
    A.set(0, 0, 2.0); A.set(0, 1, 1.0); A.set(0, 2, 1.0);
    A.set(1, 0, 1.0); A.set(1, 1, 3.0); A.set(1, 2, 2.0);
    A.set(2, 0, 1.0); A.set(2, 1, 0.0); A.set(2, 2, 0.0);
    
    std::vector<double> x_known = {1.0, 2.0, 3.0};

    std::vector<double> b(3, 0.0);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            b[i] += A.get(i, j) * x_known[j];
        }
    }
    
    std::vector<double> x = solve_slay(A, b);
    
    std::cout << "Expected: ";
    for (double val : x_known) std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "Got: ";
    for (double val : x) std::cout << val << " ";
    std::cout << std::endl;
    return 0;
}