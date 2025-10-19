#pragma once

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

constexpr double kEps = 1e-14;
constexpr size_t kWidth = 15;

class Matrix;

Matrix Identity(size_t dim);
Matrix Transpose(const Matrix& mat);

class Matrix {
public:
    Matrix(size_t rows, size_t cols) {
        data_ = std::vector<std::vector<double>>(rows,
                                                 std::vector<double>(cols, 0));
    }

    Matrix(size_t dim) {
        data_ =
            std::vector<std::vector<double>>(dim, std::vector<double>(dim, 0));
    }

    Matrix(const std::vector<std::vector<double>>& data) : data_(data) {}

    Matrix(const std::vector<double>& diag) {
        Matrix mat = Matrix(diag.size());
        for (size_t i = 0; i < diag.size(); ++i) {
            mat[i, i] = diag[i];
        }
        *this = mat;
    }

    size_t Rows() const { return data_.size(); }

    size_t Columns() const { return data_.front().size(); }

    bool operator==(const Matrix& other) const = default;

    double& operator[](size_t row, size_t col) { return data_[row][col]; }

    double operator[](size_t row, size_t col) const { return data_[row][col]; }

    std::vector<double>& operator[](size_t row) { return data_[row]; }

    Matrix operator*(const Matrix& other) const {
        Matrix res = Matrix(this->Rows(), other.Columns());
        for (size_t i = 0; i < res.Rows(); ++i) {
            for (size_t j = 0; j < res.Columns(); ++j) {
                double elem = 0.;
                for (size_t k = 0; k < this->Columns(); ++k) {
                    elem += (*this)[i, k] * other[k, j];
                }
                res[i, j] = elem;
            }
        }
        return res;
    }

    Matrix& operator*=(const Matrix& other) {
        Matrix mul = (*this) * other;
        *this = mul;
        return *this;
    }

    std::vector<double> GetDiagonal() const {
        std::vector<double> res(data_.size());
        for (size_t i = 0; i < data_.size(); ++i) {
            res[i] = data_[i][i];
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& mat) {
        for (const auto& row : mat.data_) {
            for (const auto& elem : row) {
                out << std::setw(kWidth) << elem << " ";
            }
            out << "\n";
        }
        return out;
    }

    void MultiplyAddRow(size_t bear_row, double val, size_t row) {
        for (size_t j = 0; j < this->Columns(); ++j) {
            data_[row][j] += data_[bear_row][j] * val;
        }
    }

    double FindMaxElemInColumn(size_t row, size_t col) {
        size_t idx = row;
        for (size_t i = row; i < this->Rows(); ++i) {
            if (data_[i][col] > data_[idx][col]) {
                idx = i;
            }
        }
        std::swap(data_[idx], data_[row]);
        if (std::abs(data_[row][col]) < kEps) {
            throw std::runtime_error("determinant of matrix = 0");
        }
        return data_[row][col];
    }

    Matrix FindLUDecomp() {
        Matrix res = Identity(this->Rows());
        for (size_t j = 0; j < this->Columns(); ++j) {
            double op_elem = data_[j][j];
            if (std::abs(op_elem) < kEps) {
                op_elem = this->FindMaxElemInColumn(j, j);
            }
            for (size_t i = j + 1; i < this->Rows(); ++i) {
                double rev = -data_[i][j] / op_elem;
                this->MultiplyAddRow(j, rev, i);
                res[i, j] = -rev;
            }
        }
        return res;
    }

private:
    std::vector<std::vector<double>> data_;
};

Matrix Transpose(const Matrix& mat) {
    Matrix res = Matrix(mat.Columns(), mat.Rows());
    for (size_t i = 0; i < mat.Rows(); ++i) {
        for (size_t j = 0; j < mat.Columns(); ++j) {
            res[j, i] = mat[i, j];
        }
    }
    return res;
}

Matrix Identity(size_t dim) {
    Matrix res = Matrix(dim);
    for (size_t i = 0; i < res.Rows(); ++i) {
        res[i, i] = 1;
    }
    return res;
}
