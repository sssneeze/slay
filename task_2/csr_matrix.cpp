#include <iostream>
#include <vector>

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