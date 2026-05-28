#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static constexpr double kEps = 1e-10;
static constexpr double kPi = 3.14159265358979;

std::vector<double> x, c, d, e, f, dx, res;

void Progonka(int n) {
    const int kM = n - 1;

    for (int k = 2; k <= kM; ++k) {
        double factor = c[k] / d[k - 1];
        d[k] = d[k] - e[k - 1] * factor;
        f[k] = f[k] - f[k - 1] * factor;
    }

    dx[kM] = f[kM] / d[kM];
    for (int k = kM - 1; k >= 1; --k) {
        dx[k] = (f[k] - e[k] * dx[k + 1]) / d[k];
    }

    dx[0] = 0.;
    dx[n] = 0.;
}

void Newton(int n) {
    const double kH = 2. * kPi / n;

    x.resize(n + 1);
    c.resize(n + 1);
    d.resize(n + 1);
    e.resize(n + 1);
    f.resize(n + 1);
    dx.resize(n + 1, 0.);

    res.clear();

    for (int i = 0; i <= n; ++i) {
        x[i] = 10. * (1. - i * 1. / n);
    }

    while (true) {
        double max_res = 0.0;

        for (int i = 1; i < n; ++i) {
            double fi = x[i - 1] - 2 * x[i] + x[i + 1] - 2 * kH * kH +
                        (kH / 2.) * x[i] * (x[i + 1] - x[i - 1]);

            max_res = std::max(max_res, std::abs(fi));

            f[i] = -fi;
            c[i] = 1. - (kH / 2.) * x[i];
            d[i] = -2. + (kH / 2.) * (x[i + 1] - x[i - 1]);
            e[i] = 1. + (kH / 2.) * x[i];
        }

        res.push_back(max_res);

        if (max_res < kEps) {
            break;
        }

        Progonka(n);

        for (int i = 1; i < n; ++i) {
            x[i] += dx[i];
        }
    }

    std::ofstream fout_conv("conv" + std::to_string(n) + ".txt");
    fout_conv << std::setprecision(10);
    for (double r : res) {
        fout_conv << r << "\n";
    }

    std::ofstream fout_sol("sol" + std::to_string(n) + ".txt");
    fout_sol << std::setprecision(10);
    for (int i = 0; i <= n; ++i) {
        fout_sol << i << " " << x[i] << "\n";
    }
}

int main() {
    Newton(10);
    Newton(100);
    Newton(1000);
    return 0;
}
