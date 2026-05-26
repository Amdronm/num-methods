#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

static constexpr size_t kDim = 2;

using Vec = std::array<double, kDim>;

constexpr Vec operator+(Vec a, Vec b) {
    Vec res;
    for (size_t i = 0; i < kDim; ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}
constexpr Vec operator-(Vec a, Vec b) {
    Vec res;
    for (size_t i = 0; i < kDim; ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}
constexpr Vec operator*(double a, Vec b) {
    Vec res;
    for (size_t i = 0; i < kDim; ++i) {
        res[i] = a * b[i];
    }
    return res;
}

static constexpr double kA = 0.;
static constexpr double kB = 10.;
static constexpr Vec kY0 = {2., 1.};

// change
static constexpr size_t kN = 1e8;
static constexpr Vec kY1 = {2.0000002,1.};

static constexpr double kH = (kB - kA) / kN;
std::string kPath = "output" + (kN == 200 ? "" : std::to_string(kN)) + ".csv";

// specific task
constexpr Vec F(Vec y) {
    Vec res{0};
    auto [u, v] = y;
    res[0] = 3 * u - 2 * u * v;
    res[1] = u * v - 2 * v;
    return res;
}

// Adams method
constexpr Vec Yi(Vec yi_1, Vec fi_1, Vec fi_2) {
    return yi_1 + kH * (1.5 * fi_1 - 0.5 * fi_2);
}

int main() {
    std::ofstream fout(kPath);

    std::vector<std::array<double, kDim>> ys = {kY0, kY1};
    ys.resize(kN + 1);

    Vec f_prev = F(kY0);
    for (size_t i = 2; i < ys.size(); ++i) {
        Vec f_curr = F(ys[i - 1]);
        ys[i] = Yi(ys[i - 1], f_curr, f_prev);
        f_prev = f_curr;
    }

    double x = kA;
    fout << std::setprecision(14);
    // for (auto [u, v] : ys) {
    //     fout << x << "," << u << "," << v << "\n";
    //     x += kH;
    // }
    auto [u, v] = ys.back();
    fout << u << "," << v;

    return 0;
}
