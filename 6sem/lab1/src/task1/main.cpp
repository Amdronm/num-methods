#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>

static constexpr double kEps = 1.e-12;
static constexpr size_t kIters = 1e3;

const double kRevLog16 = 1. / std::log(16.);
const double kLog16 = std::log(16.);

double F(double x) { return std::pow(16., -x) + std::log(x) * kRevLog16; }
double Df(double x) { return -std::pow(16., -x) * kLog16 + 1. / x * kRevLog16; }

const std::vector<std::pair<double, double>> kIntervals = {
    {.2182736, .3}, {.3, .45}, {.45, .54789798}};

template <typename T>
void MyPrint(const std::vector<T>& vec, std::ostream& out = std::cerr) {
    out << (vec.empty() ? -1 : static_cast<int>(vec.size())) << "\n";
    for (const T& elem : vec) {
        out << elem << " ";
    }
    out << "\n";
}

double BisectionMethod(double a, double b, std::vector<double>& pts) {
    const double kFa = F(a);
    double x = a;
    double fx = F(a);

    size_t cnt{};
    while (std::abs(fx) > kEps) {
        x = (a + b) / 2.;
        fx = F(x);

        pts.push_back(std::abs(fx));
        if (fx * kFa < 0) {
            b = x;
        } else {
            a = x;
        }

        if (cnt++ > kIters || x < 0) {
            std::cerr << "[" << a << ", " << b << "] is broken\n";
            pts.clear();
            break;
        }
    }
    return x;
}

double NewtonConstDfMethod(double a, double b, std::vector<double>& pts) {
    const double kRevMu = 1. / Df(a);
    std::cerr << "mu=" << 1. / kRevMu << "\n";

    double x = a;
    double fx = F(a);

    size_t cnt{};
    while (std::abs(fx) > kEps) {
        x -= fx * kRevMu;

        fx = F(x);
        pts.push_back(std::abs(fx));
        if (cnt++ > kIters || x < 0) {
            std::cerr << "x = " << x << ", [" << a << ", " << b
                      << "] is broken\n";
            pts.clear();
            break;
        }
        if (cnt == 1) {
            std::cerr << "x1 = " << x << "\n";
        }
    }
    return x;
}

void ProcMethod(
    std::function<double(double, double, std::vector<double>&)> method,
    std::ostream& out) {
    out << std::setprecision(12);

    std::vector<double> pts;
    for (auto [a, b] : kIntervals) {
        out << method(a, b, pts) << "\n";
        MyPrint(pts, out);
        pts.clear();
    }
}

int main() {
    std::ofstream foutb("outputB.txt");
    std::ofstream foutn("outputN.txt");

    ProcMethod(BisectionMethod, foutb);
    ProcMethod(NewtonConstDfMethod, foutn);

    return 0;
}
