#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <cmath>
#include <climits>

#define _MN_ERR 1e-12

struct range {
    size_t l, r, len;
    range(size_t l, size_t r) : l(l), r(r), len(r-l+1) {
        assert(0<=l && l<=r);
    }
};

template <typename T>
class Matrix {

public:
    // constructor with memory allocation
    Matrix(size_t r, size_t c, T x=0) : R(r), C(c) {
        M.assign(R, std::vector<T>(C, x));
    }

    // constructor from initializer_list
    Matrix(const std::initializer_list<std::vector<T>> m) : M(m) {
        R = M.size(), C = R > 0 ? M[0].size() : 0;
    }

    //  access size (rows, columns)
    bool empty() const {
        return R*C == 0;
    }

    //  access size (rows, columns)
    std::pair<size_t, size_t> size() {
        return std::make_pair(R, C);
    }

    // Algebric Operations:
    // 1 Matrix(A), 1 Scalar(b): A+I*b, A-I*b, A*b
    // 2 Matrices(A, B): A+B, A-B, A*B
    Matrix operator+=(const Matrix &rhs) {
        assert(R == rhs.R && C == rhs.C);
        for (size_t r = 0; r < R; r++) {
            for (size_t c = 0; c < C; c++) {
                M[r][c] += rhs(r, c);
            }
        }
        return *this;
    }
    Matrix operator-=(const Matrix &rhs) {
        assert(R == rhs.R && C == rhs.C);
        for (size_t r=0; r<R; r++) {
            for (size_t c=0; c<C; c++) {
                M[r][c] -= rhs(r, c);
            }
        }
        return *this;
    }
    Matrix operator*=(const Matrix &rhs) {
        assert(C == rhs.R);
        Matrix res(R, rhs.C, 0);
        for (size_t r=0; r<res.R; r++) {
            for (size_t c=0; c<res.C; c++) {
                for (size_t i=0; i<C; i++) {
                    res(r, c) += M[r][i] * rhs(i, c);
                }
            }
        }
        return res;
    }
    Matrix operator+(const Matrix &rhs) const {
        return Matrix(*this) += rhs;
    }
    Matrix operator-(const Matrix &rhs) const {
        return Matrix(*this) -= rhs;
    }
    Matrix operator*(const Matrix &rhs) const {
        return Matrix(*this) *= rhs;
    }
    Matrix operator*=(const T &rhs) {
        for(size_t r=0; r<R; r++) {
            for(size_t c=0; c<C; c++) {
                M[r][c] *= rhs;
            }
        }
        return *this;
    }
    Matrix operator*(const T &rhs) const {
        return Matrix(*this) *= rhs;
    }
    Matrix operator-() const {
        return Matrix(R, C, 0) - *this;
    }

    // L2 norm of A-B (ie ||A-B||^2) must be lesser than _MN_ERR (=1e-12)
    bool operator==(const Matrix &rhs) const {
        if (empty() || R != rhs.R || C != rhs.C) {
            return false;
        }
        double L2 = 0;
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < C; j++) {
                L2 += (M[i, j] - rhs(i, j)) * (M[i, j] - rhs(i, j));
            }
        }
        L2 /= R*C;
        return L2 < _MN_ERR;
    }
    bool operator!=(const Matrix &rhs) const {
        return !(*this == rhs);
    }

    // Element Numbering (0 based, Column  Major):
    //    C 0  1  2   3
    //  R
    //  0   0  3  6   9
    //  1   1  4  7  10
    //  2   2  5  8  11

    // access data operators
    T& operator[](size_t i) {
        return M[i % R][i / R];
    }
    T& operator()(size_t r, size_t c) {
        assert(0 <= r && r < R && 0 <= c && c < C);
        return M[r][c];
    }
    T operator()(size_t r, size_t c) const {
        assert(0 <= r && r < R && 0 <= c && c < C);
        return M[r][c];
    }
    Matrix operator()(const range &row_range, const range &col_range) const {
        assert(row_range.r < R && col_range.r < C);
        Matrix<T> res(row_range.len, col_range.len, 0);
        for (size_t r = 0; r < row_range.len; r++) {
            for (size_t c = 0; c < col_range.len; c++)
                res(r, c) = M[row_range.l + r][col_range.l + c];
        }
        return res;
    }

    // compute transpose
    Matrix t() {
        Matrix res(C, R);
        for(size_t r=0; r<R; r++) {
            for(size_t c=0; c<C; c++) {
                res(c, r) += M[r][c];
            }
        }
        return res;
    }

    // compute dth minor
    Matrix compute_minor(int d) {
        Matrix res(R, C);
        for (int i = 0; i < d; i++) {
            res(i, i) = 1;
        }
        for (int r = d; r < R; r++) {
            for (int c = d; c < C; c++) {
                res(r, c) = M[r][c];
            }
        }
        return res;
    }

    // return c-th column of m
    std::vector<T> column(int c) const {
        std::vector<T> v(R);
        for (size_t r = 0; r < R; r++) {
            v[r] = M[r][c];
        }
        return v;
    }

    // return r-th row of m
    std::vector<T> row(int r) const {
        return M[r];
    }

    // Matrix m to exponent n (Binary Exponentiation, O(logN))
    template <typename U>
    friend Matrix<U> pow(Matrix<U> m, int n);

    // To permit: 5*M, without this only M*5 is permitted
    template <typename U>
    friend Matrix<U> operator*(const U &lhs, const Matrix<U> &rhs);

private:
    size_t R, C;
    std::vector<std::vector<T>> M;
};

template <typename T>
Matrix<T> operator*(const T &lhs, const Matrix<T> &rhs) {
    return rhs * lhs;
}

template <typename T>
std::ostream& operator<<(std::ostream &cout, Matrix<T> m) {
    size_t R = m.size().first, C = m.size().second;
    for(size_t r = 0; r < R; r++) {
        std::cout << "[";
        for(size_t c = 0; c < C; c++) {
            cout << m(r, c);
            if (c+1 != C) {
                cout << " ";
            }
            else {
                cout << "]" << std::endl;
            }
        }
    }
    return cout;
}

// Zero/null Matrix of size r*c
template <typename T>
Matrix<T> zeros(size_t r, size_t c) {
    return Matrix<T>(r, c, 0);
}

// Identity Matrix of size r*r
template<typename T>
Matrix<T> eye(size_t r) {
    Matrix<T> res(r, r, 0);
    for(size_t i=0; i<r; i++) {
        res(i, i) = 1;
    }
    return res;
}

template<typename T>
Matrix<T> pow(Matrix<T> m, int n) {
    assert(m.R == m.C && n >= 0);
    if(n == 0) {
        return eye<T>(m.R);
    }
    Matrix<T> res = m;
    for(int i=n; i>0; i>>=1) {
        res = (i&1) ? res*m : res;
        m *= m;
    }
    return res;
}

// Norm of Matrix m
// Options available:
// L1: Manhatten Norm
// L2: Eucledian Norm (under development)
// LF: Frobenius Norm
// Linf: Maximum Norm
template<typename T>
double norm(Matrix<T> m, std::string NORM = "LF") {
    double res = 0;
    size_t R = m.size().first, C = m.size().second;
    if (NORM == "L1") {
        for(size_t c = 0; c < C; c++) {
            double cur = 0;
            for(size_t r = 0; r < R; r++) {
                cur += abs(m(r, c));
            }
            res = std::max(res, cur);
        }
    }
    else if (NORM == "LF"){
        for(size_t r = 0; r < R; r++) {
            for(size_t c = 0; c < C; c++) {
                res += m(r, c) * m(r, c);
            }
        }
        res = sqrt(res);
    }
    else if (NORM == "Linf") {
        for(size_t r = 0; r < R; r++) {
            double cur = 0;
            for(size_t c = 0; c < C; c++) {
                cur += abs(m(r, c));
            }
            res = std::max(res, cur);
        }
    }
    else {
        assert(false);
    }
    return res;
}

// kth norm of vector v
// k != inf: pow(sum(v^k), 1/k)
// k == inf: max element in v
template<typename T>
double norm(std::vector<T> v, std::string NORM = "L2") {
    double res = 0;
    size_t n = v.size();
    if(NORM == "L1") {
        for(size_t i = 0; i < n; i++) {
            res += std::abs(v[i]);
        }
    }
    else if(NORM == "L2") {
        for(size_t i = 0; i < n; i++) {
            res += v[i] * v[i];
        }
        res = std::sqrt(res);
    }
    else if (NORM == "Linf") {
        for(size_t i = 0; i < n; i++) {
            res = std::max(res, v[i]);
        }
    }
    else {
        assert(false);
    }
    return res;
}

// MaxVal - MinVal over axis
// axis 1: range of each column
// axis 2: range of each row
template<typename T>
std::vector<T> range(Matrix<T> m, int axis = 1) {
    assert(axis == 1 || axis == 2);
    size_t R = m.size().first, C = m.size().second;
    std::vector<T> res(axis == 1 ? C : R);

    if(axis == 1) {
        for(size_t c = 0; c < C; c++) {
            T mnVal = std::numeric_limits<T>::max();
            T mxVal = std::numeric_limits<T>::lowest();
            for(size_t r = 0; r < R; r++) {
                mnVal = std::min(mnVal, m(r, c));
                mxVal = std::max(mxVal, m(r, c));
            }
            res[c] = mxVal - mnVal;
        }
    }
    else {
        for(size_t r = 0; r < R; r++) {
            T mnVal = std::numeric_limits<T>::max();
            T mxVal = std::numeric_limits<T>::lowest();
            for(size_t c = 0; c < C; c++) {
                mnVal = std::min(mnVal, m(r, c));
                mxVal = std::max(mxVal, m(r, c));
            }
            res[r] = mxVal - mnVal;
        }
    }
    return res;
}

// divide data by factor
template<typename T>
std::vector<T> rescale(const std::vector<T> &v, T factor=0) {
    if(factor == 0) {
        factor = norm(v);
    }
    assert(factor != 0);
    std::vector<T> res(v);
    for (int i = 0; i < v.size(); i++) {
        res[i] /= factor;
    }
    return res;
}

#endif