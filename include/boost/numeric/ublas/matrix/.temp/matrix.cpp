//****************************//   RELAXXPLS   //****************************//
#include <bits/stdc++.h>
using namespace std;
// #include <vector>
// #include <utility>
// #include <iostream>
// #include <algorithm>
// #include <cassert>
// #include <iterator>

/******************************** Debug Begin ********************************/
#define ts to_string
#ifdef RELAXXPLS
    #define cerr cout
    #define deb(x...) cerr<<#x<<": ", _deb(x)
#else
    #define deb(x...)
#endif
string ts(string c) { return string(c); }
string ts(char ch)  { return string(1, ch); }
string ts(bool b)   { return b ? "true" : "false"; }
template<size_t N> string ts(bitset<N> B) {
    string S=""; for(int i=0; i<(int)N; i++) S += (char)'0'+B[i]; return S; }
template<typename T, typename U> string ts(pair<T, U> P) {
    return "{"+ts(P.first)+","+ts(P.second)+"}"; }
template<typename T> string ts(T A) { string S="["; bool F=0;
    for(auto x: A) { if(F) S+=" "; F=1; S += ts(x); } return S+"]"; }
template<typename T, size_t N> string ts(T (&A)[N]) { string S="["; bool F=0;
    for(int i=0; i<N; i++) { if(F) S+=" "; F=1; S+=ts(A[i]); } return S+"]"; }
void _deb() { cerr<<endl; }
template<typename T, typename... U> void _deb(const T& t, const U&... u) {
    cerr<<ts(t); if(sizeof...(u)) cout<<", "; _deb(u...); }
/********************************* Debug Eng *********************************/

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

//============================//   Code End   //=============================//
int main() {
    matrix<int> M(2, 9, 1);
    M[7] = 7;
    M(0, 1) = 2;
    std::cout << M << std::endl;
    matrix<int> M2 = M(range(1, 1), range(2, 6));
    std::cout << M2 << std::endl;

    return 0;
}

// #pragma once
// #include <vector>
// #include <utility>
// #include <iostream>
// #include <algorithm>
// #include <cassert>
// #include <iterator>
// #include "matrix.hpp"

// namespace boost { namespace numeric { namespace ublas {

// template<typename T>
// class matrix {
// public:
// 	typedef int size_type;
// 	typedef T value_type;

// 	matrix(size_type r, size_type c, T x) : R(r), C(c) {
// 		M.assign(R, std::vector<T>(C, x));
// 	}
// 	~matrix() {}

// 	bool empty() {
// 		return R*C == 0;
// 	}

// 	matrix operator+=(const matrix &rhs) {
// 		assert(R==rhs.R && C==rhs.C);
//         for(size_type r=0; r<R; r++) {
//         	for(size_type c=0; c<C; c++) {
//         		M[r][c] += rhs.M[r][c];
//         	}
//         }
//         return *this;
// 	}
//     matrix operator-=(const matrix &rhs) {
//     	assert(R==rhs.R && C==rhs.C);
//         for(size_type r=0; r<R; r++) {
//         	for(size_type c=0; c<C; c++) {
//         		M[r][c] -= rhs.M[r][c];
//         	}
//         }
//         return *this;
//     }
//     matrix operator*=(const matrix &rhs) {
//     	assert(C == rhs.R);
//         matrix res(R, rhs.C, 0);
//         for(size_type r=0; r<res.R; r++) {
//         	for(size_type c=0; c<res.C; c++) {
//             	for(size_type i=0; i<C; i++) {
//             		res.M[r][c] += M[r][i] * rhs.M[i][c];
//             	}
//         	}
//         }
//         return res;
//     }
//     matrix operator*=(const T &rhs) {
//         for(size_type r=0; r<R; r++) {
//         	for(size_type c=0; c<C; c++) {
//         		M[r][c] *= rhs;
//         	}
//         }
//         return *this;
//     }
//     matrix operator-() {
//     	return matrix(R, C, 0) - *this;
//     }
//     matrix operator+() {
//     	return *this;
//     }

//     friend matrix operator+(matrix a, const matrix &b) {
//         return a += b;
//     }
//     friend matrix operator-(matrix a, const matrix &b) {
//         return a -= b;
// 	}
//     friend matrix operator*(matrix a, const matrix &b) {
//         return a *= b;
//     }

//     // friend matrix pow(matrix m, int n) {
//     //     matrix res = eye(m.R);
//     //     for(int i=n; i>0; i>>=1) {
//     //     	res = (i&1) ? res*m : res;
//     //     	m *= m;
//     //     }
//     //     return res; }
// private:
// 	size_type R, C;
// 	std::vector<std::vector<value_type>> M;
// };

// // Zero/null matrix of size r*c
// template<typename T>
// matrix<T> zeros(size_t r, size_t c) {
// 	return matrix<T>(r, c, 0);
// }

// // Identity matrix of size r*r
// template<typename T>
// matrix<T> eye(size_t r) {
// 	matrix<T> res(r, r, 0);
// 	for(int i=0; i<r; i++) {
// 		res(i, i) = 1;
// 	}
// 	return res;
// }

// }}} // Namespace end