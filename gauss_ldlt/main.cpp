#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <vector>

#include "matrix.h"

Matrix ReadMatrix(std::istream& fin) {
    if (fin.bad()) {
        std::cerr << "Wrong file" << std::endl;
        return 0;
    }
    size_t n = 0;
    fin >> n;

    Matrix mat_a = Matrix(n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            fin >> mat_a[i, j];
        }
    }
    return mat_a;
}

// solves mat*x=vec system with bot triangle mat
std::vector<double> SolveBotTriangle(const Matrix& mat,
                                     const std::vector<double>& vec) {
    std::vector<double> sol;
    sol.reserve(vec.size());
    sol.push_back(vec.front());
    for (size_t i = 1; i < vec.size(); ++i) {
        double x = vec[i];
        for (size_t j = 0; j < i; ++j) {
            x -= mat[i, j] * sol[j];
        }
        sol.push_back(x);
    }
    return sol;
}

// solves mat*x=vec system with top triangle mat
std::vector<double> SolveTopTriangle(const Matrix& mat,
                                     const std::vector<double>& vec) {
    std::vector<double> sol(vec.size());
    sol.back() = vec.back();
    for (int i = vec.size() - 2; i >= 0; --i) {
        double x = vec[i];
        for (size_t j = i; j < vec.size(); ++j) {
            x -= mat[i, j] * sol[j];
        }
        sol[i] = x;
    }
    return sol;
}

void MulVector(std::vector<double>& vec, const std::vector<double>& vec_by) {
    for (size_t i = 0; i < vec.size(); ++i) {
        vec[i] *= vec_by[i];
    }
}

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& other) {
    for (const auto& elem : other) {
        out << std::setw(kWidth) << elem << " ";
    }
    out << "\n";
    return out;
}

// returns matrix maximum norm
double FindMaxNorm(const Matrix& mat) {
    double y = 0.;
    for (size_t i = 0; i < mat.Rows(); ++i) {
        double norm = 0.;
        for (size_t j = 0; j < mat.Columns(); ++j) {
            norm += std::fabs(mat[i, j]);
        }
        y = std::max(norm, y);
    }
    return y;
}

void SolveFast(std::istream& in, std::ostream& out);

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);

    std::ifstream fin;
    std::ofstream fout;
    std::istream* in_ptr;
    std::ostream* out_ptr;

    if (argc != 3) {
        fin.open("input.txt");
        in_ptr = &fin;
        out_ptr = &std::cout;
    } else {
        fin.open(argv[1]);
        in_ptr = &fin;

        fout.open(argv[2]);
        out_ptr = &fout;

        if (!fin.is_open() || !fout.is_open()) {
            std::cerr << "wrong filename(s)" << std::endl;
            return 1;
        }
    }
    auto& in = *in_ptr;
    auto& out = *out_ptr;

    SolveFast(in, out);

    return 0;
}

void SolveFast(std::istream& in, std::ostream& out) {
    Matrix mat_a = ReadMatrix(in);
    double norm_a = FindMaxNorm(mat_a);

    Matrix mat_l = mat_a.FindLUDecomp();
    auto diag_d = mat_a.GetDiagonal();

    out << "Matrix L :\n" << mat_l << "\nDiagonal of Matrix D : \n" << diag_d;
    auto diag_d_rev = diag_d;
    for (auto& elem : diag_d_rev) {
        elem = 1. / elem;
    }

    /// here we find A^-1 solving
    /// n eqations Ax=e_i i = 1,2,...,n
    /// where e_i = (0, 0, ..., 1, ... 0) 1 in i position
    /// first Ly_i=e_i
    /// y *= D^-1
    /// L^t*x_i=y_i
    /// x_i is i column of A^-1
    /// thus we solved A*A^-1=E column by column
    size_t dim = diag_d.size();
    Matrix rev_a(dim);
    auto mat_lt = Transpose(mat_l);
    for (size_t i = 0; i < dim; ++i) {
        std::vector<double> vec_e(dim, 0.);
        vec_e[i] = 1.;
        auto y = SolveBotTriangle(mat_l, vec_e);
        MulVector(y, diag_d_rev);
        y = SolveTopTriangle(mat_lt, y);
        rev_a[i] = y;
    }

    out << "\nReversed Matrix A : \n" << rev_a;

    out << "\nCondition number = "
        // << std::setprecision(15)
        << norm_a * FindMaxNorm(rev_a) << std::endl;
}
