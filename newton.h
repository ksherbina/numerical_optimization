#ifndef NEWTON_H
#define NEWTON_H
 
#include <valarray>
using std::valarray;

#include <string>
using std::string;

valarray<double> newton(valarray<double> x0, double (*func)(valarray<double>), valarray<double> (*grad)(valarray<double>),
                                valarray<double> (*hess)(valarray<double>), double epsilon, string newton_type, string linsearch_type, int MAX_ITER);
 
#endif
