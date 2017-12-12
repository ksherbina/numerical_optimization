#ifndef DAMPED_BFGS_H
#define DAMPED_BFGS_H
 
#include <valarray>

valarray<double> damped_bfgs(valarray<double> B0, valarray<double> delta, valarray<double> gamma);
 
#endif
