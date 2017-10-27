#ifndef NEWTON_DIRECTION_H
#define NEWTON_DIRECTION_H
 
#include <valarray>
#include "cholesky.h"
#include "solve_linear_system.h"
using std::valarray;

valarray<double> newton_direction(double (*f)(valarray<double>),valarray<double> (*gradient)(valarray<double>),
  valarray<double> (*hessian)(valarray<double>),valarray<double> x0, int dim, double epsilon, int modify_diagonal, int linesearch);
 
#endif
