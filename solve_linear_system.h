#ifndef SOLVE_LINEAR_SYSTEM_H
#define SOLVE_LINEAR_SYSTEM_H
 
#include <valarray>
using std::valarray;

valarray<double> solve_linear_system(valarray<double> L, valarray<double> D, valarray<double> b, int n);
 
#endif
