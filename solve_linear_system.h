#ifndef SOLVE_LINEAR_SYSTEM_H
#define SOLVE_LINEAR_SYSTEM_H
 
#include <valarray>
#include "cholesky.h"
using std::valarray;

valarray<double> solve_linear_system(CholeskyFactors factors, valarray<double> b, int n);
 
#endif
