#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// ================= ЗАДАЧА 1 =================

double F1(double x) { return pow(16, -x) + log(x) / log(16.0); }

double Df1(double x) {
    return -pow(16, -x) * log(16.0) + 1.0 / (x * log(16.0));
}

void SolveTask1(double a, double b, int root_idx) {
    double eps = 1e-12;
    ofstream out("task1_root" + to_string(root_idx) + ".csv");
    out << "Iter,Bisection_Res,Newton_Res\n";

    // Метод бисекции
    double left = a;
    double right = b;
    double c;
    vector<double> err_bisect;
    while (true) {
        c = (left + right) / 2.0;
        double fc = F1(c);
        err_bisect.push_back(abs(fc));
        if (abs(fc) < eps) break;
        if (F1(left) * fc < 0)
            right = c;
        else
            left = c;
    }

    // Метод Ньютона с постоянной производной (x0 берем как середину отрезка)
    double x0 = (a + b) / 2.0;
    double const_df = Df1(x0);
    double xk = x0;
    vector<double> err_newton;
    while (true) {
        double fk = F1(xk);
        err_newton.push_back(abs(fk));
        if (abs(fk) < eps) break;
        xk = xk - fk / const_df;
    }

    // Запись в файл
    size_t max_iter = max(err_bisect.size(), err_newton.size());
    for (size_t i = 0; i < max_iter; ++i) {
        out << i << ",";
        if (i < err_bisect.size())
            out << scientific << err_bisect[i] << ",";
        else
            out << ",";
        if (i < err_newton.size()) out << scientific << err_newton[i];
        out << "\n";
    }
    out.close();
    cout << "Root " << root_idx << " found: " << setprecision(15) << c
         << " (Bisect) / " << xk << " (Newton)\n";
}

// ================= ЗАДАЧА 2 =================

void SolveTask2(int n) {
    double h = 2 * M_PI / n;
    vector<double> x(n + 1);

    // Начальное приближение
    for (int i = 0; i <= n; ++i) {
        x[i] = 10.0 * (1.0 - (double)i / n);
    }

    double eps = 1e-10;
    vector<double> A(n);
    vector<double> B(n);
    vector<double> C(n);
    vector<double> F(n);
    vector<double> alpha(n);
    vector<double> beta(n);
    vector<double> dx(n + 1, 0.0);

    ofstream out_conv("task2_conv_" + to_string(n) + ".csv");
    out_conv << "Iter,NormF\n";

    int iter = 0;
    while (true) {
        double normF = 0.0;

        // Формирование СЛАУ для метода прогонки
        for (int i = 1; i < n; ++i) {
            A[i] = 1.0 - (h / 2.0) * x[i];
            C[i] = -2.0 + (h / 2.0) * (x[i + 1] - x[i - 1]);
            B[i] = 1.0 + (h / 2.0) * x[i];

            double fi = x[i - 1] - 2 * x[i] + x[i + 1] - 2 * h * h +
                        (h / 2.0) * x[i] * (x[i + 1] - x[i - 1]);
            F[i] = -fi;

            normF = max(
                normF, abs(fi));  // Чебышевская (векторная бесконечность) норма
        }

        out_conv << iter << "," << scientific << normF << "\n";

        if (normF < eps) {
            break;
        }

        // Прямой ход прогонки
        alpha[1] = -B[1] / C[1];
        beta[1] = F[1] / C[1];
        for (int i = 2; i < n; ++i) {
            double denom = A[i] * alpha[i - 1] + C[i];
            alpha[i] = -B[i] / denom;
            beta[i] = (F[i] - A[i] * beta[i - 1]) / denom;
        }

        // Обратный ход прогонки
        dx[n - 1] = beta[n - 1];
        for (int i = n - 2; i >= 1; --i) {
            dx[i] = alpha[i] * dx[i + 1] + beta[i];
        }

        // Обновление x
        for (int i = 1; i < n; ++i) {
            x[i] += dx[i];
        }
        iter++;
    }
    out_conv.close();

    ofstream out_sol("task2_sol_" + to_string(n) + ".csv");
    out_sol << "i,xi\n";
    for (int i = 0; i <= n; ++i) {
        out_sol << i << "," << x[i] << "\n";
    }
    out_sol.close();
}

int main() {
    cout << "--- Task 1 ---" << endl;
    SolveTask1(0.1, 0.3, 1);
    SolveTask1(0.3, 0.4, 2);
    SolveTask1(0.4, 0.6, 3);

    cout << "\n--- Task 2 ---" << endl;
    SolveTask2(10);
    SolveTask2(100);
    SolveTask2(1000);

    cout << "\nCalculations done. Run the Python script to plot charts."
         << endl;
    return 0;
}
