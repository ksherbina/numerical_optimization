#ifndef CHOLESKY_RANKTWO_H
#define CHOLESKY_RANKTWO_H
 
#include <valarray>
#include "cholesky.h"
using std::valarray;

CholeskyFactors cholesky_ranktwo(valarray<double> v, valarray<double>u, CholeskyFactors chol0);
 
#endif
