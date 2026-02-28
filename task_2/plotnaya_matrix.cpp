#include <iostream>
#include <vector>

class PlotMatrix {
public:
    PlotMatrix(size_t m, size_t n) : m_(m), n_(n) {
        data.resize(m * n, 0.0);
    }

    double get(size_t i, size_t j) const {
        return data[i * n_ + j];
    }

    double set(size_t i, size_t j, size_t value) {
        data[i * n_ + j] = value;
    }

    void print() const {
        for (size_t i = 0; i < n_; i++) {
            for (size_t j = 0; j < m_; j++) {
                std::cout << data[i * n_ + j] << ' ';
            }
            std::cout << std::endl;
        }        
    }

    //умножение на вектор
    std::vector<double> mul_vec(const std::vector<double>& v) const {
        std::vector<double> res(m_, 0.0);

        for (size_t i = 0; i < m_; i++) {
            double sum = 0.0;
            for (size_t j = 0; j < n_; j++) {
                sum += data[i * n_ + j] * v[i];
            }
            res[i] = sum;
        }

        return res;
    }

private:
    size_t m_; //строки
    size_t n_; //столбцы
    std::vector<double> data;
};