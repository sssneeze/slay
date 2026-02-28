#include <iostream>
#include <vector>
#include <tuple>
#include <chrono>
#include <cstdlib>


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


class CSRMatrix {
public:
CSRMatrix(size_t m, size_t n, 
          const std::vector<std::tuple<int, int, double>>& dok) : m_(m), n_(n) {

    std::vector<int> count_value_rows(m, 0);
    for (const auto& elem : dok) {
        int i = std::get<0>(elem);
        double val = std::get<2>(elem);
        
        if (val != 0.0 && i >= 0 && i < m) {
            count_value_rows[i]++;
        }
    }
    rows_.resize(m + 1, 0);
    for (size_t i = 0; i < m; ++i) {
        rows_[i + 1] = rows_[i] + count_value_rows[i];
    }
    
    std::vector<int> current_pos = rows_;
    
    size_t num_of_elems = rows_[m];
    values_.resize(num_of_elems);
    cols_.resize(num_of_elems);
    
    for (const auto& elem : dok) {
        int i = std::get<0>(elem);
        int j = std::get<1>(elem);
        double val = std::get<2>(elem);

        if (val == 0.0 || i < 0 || i >= m || j < 0 || j >= n) {
            continue;
        }
        
        int pos = current_pos[i];
        values_[pos] = val;
        cols_[pos] = j;
        
        current_pos[i]++;
    }
}

    //умножение на вектор
    std::vector<double> mul_vec(const std::vector<double>& v) const {
        std::vector<double> res(m_, 0.0);

        for (size_t i = 0; i < m_; i++) {
            double sum = 0.0;

            int start = rows_[i];
            int end = rows_[i + 1];
            
            for (int j = 0; j < end; j++) {
                int col = cols_[i];
                double value = values_[j];
                sum += value * v[col];
            }
            res[i] = sum;
        }

        return res;
    }

    void print() const {
        std::cout << "Values: ";
        for (double v : values_) std::cout << v << " ";
        std::cout << std::endl;
        
        std::cout << "Cols: ";
        for (int c : cols_) std::cout << c << " ";
        std::cout << std::endl;
        
        std::cout << "Rows ptr: ";
        for (int r : rows_) std::cout << r << " ";
        std::cout << std::endl;
    }

private:
    size_t m_; //строки
    size_t n_; //столбцы
    std::vector<double> values_;
    std::vector<int> cols_;
    std::vector<int> rows_;
};


double measure_time(auto func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    return elapsed.count();
}


std::vector<std::tuple<int, int, double>> create_csr_matrix(
    size_t m, size_t n, size_t nn_per_row) 
{
    std::vector<std::tuple<int, int, double>> dok;
    for (int i = 0; i < m; i++) {
        for (int k = 0; k < nn_per_row; k++) {
            int col = rand() % n;
            double value = static_cast<double>(rand() % 100) / 10.0 + 1.0;
            dok.push_back({i, col, value});
        }
    }
    return dok;
}



int main() {
    srand(42);
    const size_t N = 100;

    std::vector<size_t> not_null_num = {5, 10, 50};

    for (size_t nn_per_row : not_null_num) {
        std::vector<std::tuple<int, int, double>> dok_matrix = 
            create_csr_matrix(N, N, nn_per_row);

        
        CSRMatrix csr_matrix(N, N, dok_matrix);

        PlotMatrix plot_matrix(N, N);
        for (const auto& elem : dok_matrix) {
            int i = std::get<0>(elem);
            int j = std::get<1>(elem);
            double val = std::get<2>(elem);
            plot_matrix.set(i, j, val);
        }

        std::vector<double> test_vec(N, 1.0);

        
        double plot_time = 0.0;
        plot_time += measure_time([&]() {
            plot_matrix.mul_vec(test_vec);
        });
        plot_time = plot_time / 3.0;

        
        double csr_time = 0.0;
        csr_time += measure_time([&]() {
            csr_matrix.mul_vec(test_vec);
        });
        csr_time = csr_time/ 3.0;

        std::cout << "Ненулевых в строке: " << nn_per_row << std::endl;
        std::cout << "Время у плотной: " << plot_time << " мс" << std::endl;
        std::cout << "Время у разреженной: " << csr_time << " мс" << std::endl;
    }

    return 0;
}