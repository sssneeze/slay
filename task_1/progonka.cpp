#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> progonka(vector<double> main_d, vector<double> lower_d, vector<double> upper_d,
                        vector<double> column, int N) {
        if (main_d.size() != N || column.size() != N ||
        (N > 1 && (lower_d.size() != N - 1 || upper_d.size() != (size_t)N - 1)))  {
            throw runtime_error("Неправильные размеры");
        }

        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            if (i > 0) sum += fabs(lower_d[i - 1]);
            if (i < N - 1) sum += fabs(upper_d[i]);
            if (fabs(main_d[i]) < sum) {
                throw runtime_error("Не выплнено условие нестрогого диагонального преобладания");
            }
        }

        vector<double> alpha(N-1);
        vector<double> beta(N);
        vector<double> X(N);
        alpha[0] = - upper_d[0] / main_d[0];
        beta[0] = column[0] / main_d[0];
        for (int i = 1; i < N - 1; i++) {
            alpha[i] = - upper_d[i] / (lower_d[i-1] * alpha[i-1] + main_d[i]);
            beta[i] = (column[i] - lower_d[i-1] * beta[i-1]) / (lower_d[i-1] * alpha[i-1] + main_d[i]);
        }
        X[N-1] = (column[N-1] - lower_d[N-2] * beta[N-2]) / (lower_d[N-2] * alpha[N-2] + main_d[N-1]);
        beta[N-1] = X[N-1];

        for (int i = N - 2; i >= 0; i--) {
            X[i] = alpha[i] * X[i+1] + beta[i];
        }

        return X;
}

int main() {
    vector<double> main_d = {4, 4, 4, 4, 4};
    vector<double> lower_d = {-1, -1, -1, -1};
    vector<double> upper_d = {-1, -1, -1, -1};
    vector<double> column = {2, 4, 6, 8, 16};
    //vector<double> column = {3, 2, 2, 2, 3};

    vector<double> result_theor = {1, 2, 3, 4, 5};
    //vector<double> result_theor = {1, 1, 1, 1, 1};

    int N = 5;
    vector<double> result_exp = progonka(main_d, lower_d, upper_d, column, N);
    for (int i = 0; i < N; i++) {
        cout << result_theor[i] << ' ';
    }
    cout << endl;
    for (int i = 0; i < N; i++) {
        cout << result_exp[i] << ' ';
    }
    cout << endl;
}