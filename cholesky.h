#ifndef CHOLESKY_H
#define CHOLESKY_H
 
#include <valarray>
using std::valarray;

struct CholeskyFactors {
    valarray<double> lower_triangular;
    valarray<double> diagonal;
};

CholeskyFactors cholesky(valarray<double> A, int n);
 
#endif
