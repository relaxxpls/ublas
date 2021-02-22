#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <iterator>

/**
 *  Operations:
 *  1 Matrix(A), 1 Scalar(b): A+I*b, A-I*b, A*b
 *  2 Matrices(A, B): A+B, A-B, A*B
 */

struct range {
    size_t l, r, len;
    range(size_t l, size_t r) : l(l), r(r), len(r-l+1) {
        assert(0<=l && l<=r);
    }
};


template <typename T>
class matrix {
public:
	matrix(size_t r, size_t c, T x) : R(r), C(c) {
        M.assign(R, std::vector<T>(C, x));
    }
	// ~matrix();

	bool empty() const {
        return R*C == 0;
    }
    std::pair<size_t, size_t> size() {
        return std::make_pair(R, C);
    }

    matrix operator+=(const matrix &rhs);
    matrix operator-=(const matrix &rhs);
    matrix operator*=(const matrix &rhs);
    matrix operator+=(const T &rhs);
    matrix operator-=(const T &rhs);
    matrix operator*=(const T &rhs);
    matrix operator+(const matrix &rhs) const;
    matrix operator-(const matrix &rhs) const;
    matrix operator*(const matrix &rhs) const;
    matrix operator+(const T &rhs) const;
    matrix operator-(const T &rhs) const;
    matrix operator*(const T &rhs) const;
    matrix operator-();
    matrix operator+();
    bool operator==(const T &rhs) const;
    bool operator!=(const T &rhs) const;


    // Element Numbering (Column  Major):
    //    C 0  1  2  3
    //  R
    //  0   0  3  6  9
    //  1   1  4  7 10
    //  2   2  5  8 11
    T& operator[](size_t i) {
        return M[i % R][i / R];
    }
    T& operator()(size_t r, size_t c) {
        return M[r][c];
    }
    matrix operator()(const range &row_range, const range &col_range) const {
        assert(row_range.r < R && col_range.r < C);
        matrix<T> res(row_range.len, col_range.len, 0);
        for (size_t r = 0; r < row_range.len; r++) {
            for (size_t c = 0; c < col_range.len; c++)
                res(r, c) = M[row_range.l + r][col_range.l + c];
        }
        return res;
    }

    template <typename TT>
    friend matrix<TT> pow(matrix<TT> m, int n);

private:
	size_t R, C;
	std::vector<std::vector<T>> M;
};

template <typename TT>
std::ostream& operator<<(std::ostream &cout, matrix<TT> m) {
    std::pair<size_t, size_t> p = m.size();
    for(size_t r = 0; r < p.first; r++) {
        std::cout << "[";
        for(size_t c = 0; c < p.second; c++) {
            cout << m(r, c);
            if (c+1 != p.second) {
                cout << " ";
            } else {
                cout << "]" << std::endl;
            }
        }
    }
    return cout;
}

// Zero/null matrix of size r*c
template <typename T>
matrix<T> zeros(size_t r, size_t c) {
	return matrix<T>(r, c, 0);
}

// Identity matrix of size r*r
template<typename T>
matrix<T> eye(size_t r) {
	matrix<T> res(r, r, 0);
	for(int i=0; i<r; i++) {
		res(i, i) = 1;
	}
	return res;
}

#endif