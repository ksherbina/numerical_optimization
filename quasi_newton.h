#ifndef QUASI_NEWTON_H
#define QUASI_NEWTON_H
 
#include <valarray>
using std::valarray;

#include <string>
using std::string;

valarray<double> quasi_newton(valarray<double> x0, double (*func)(valarray<double>), valarray<double> (*grad)(valarray<double>),
                                double epsilon, string newton_type, string linsearch_type, int MAX_ITER);
 
#endif
