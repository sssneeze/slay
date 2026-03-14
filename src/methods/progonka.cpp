#include "methods/progonka.hpp"

namespace Methods {
    std::vector<double> thomas_algorithm(const std::vector<double>& main_d,
                                         const std::vector<double>& lower_d,
                                         const std::vector<double>& upper_d,
                                         const std::vector<double>& rhs,
                                         int N) {
        // Проверка размеров
        if (main_d.size() != static_cast<size_t>(N) || 
            rhs.size() != static_cast<size_t>(N)) {
            throw std::runtime_error("Неправильные размеры main_d или rhs");
        }
        
        if (N > 1 && (lower_d.size() != static_cast<size_t>(N - 1) || 
                      upper_d.size() != static_cast<size_t>(N - 1))) {
            throw std::runtime_error("Неправильные размеры lower_d или upper_d");
        }

        // Проверка диагонального преобладания
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            if (i > 0) sum += std::fabs(lower_d[i - 1]);
            if (i < N - 1) sum += std::fabs(upper_d[i]);
            
            if (std::fabs(main_d[i]) < sum) {
                throw std::runtime_error(
                    "Не выполнено условие нестрогого диагонального преобладания");
            }
        }

        // Прямой ход
        std::vector<double> alpha(N - 1);
        std::vector<double> beta(N);
        std::vector<double> x(N);
        
        alpha[0] = -upper_d[0] / main_d[0];
        beta[0] = rhs[0] / main_d[0];
        
        for (int i = 1; i < N - 1; i++) {
            double denom = lower_d[i - 1] * alpha[i - 1] + main_d[i];
            alpha[i] = -upper_d[i] / denom;
            beta[i] = (rhs[i] - lower_d[i - 1] * beta[i - 1]) / denom;
        }
        
        // Обратный ход
        x[N - 1] = (rhs[N - 1] - lower_d[N - 2] * beta[N - 2]) / 
                   (lower_d[N - 2] * alpha[N - 2] + main_d[N - 1]);
        beta[N - 1] = x[N - 1];

        for (int i = N - 2; i >= 0; i--) {
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }

        return x;
    }
}