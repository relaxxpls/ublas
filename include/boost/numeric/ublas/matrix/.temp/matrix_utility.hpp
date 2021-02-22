#ifndef _MATRIX_UTILITY_U_
#define _MATRIX_UTILITY_H_

#include <iostream>

template <typename T>
class Size {
public:
    unsigned int r, c;
    Size(unsigned int R, unsigned int C) : r(R), c(C) {}
    unsigned int count() const {
        return r*c;
    }
    Size swap() const {
        return Size(c, r);
    }
    bool operator==(const Size &rhs) {
        return r==rhs.r && c==rhs.c;
    }
    bool operator!=(const Size &rhs) {
        return !(*this==rhs);
    }
};

class Position {

};

class Range {

};

#endif