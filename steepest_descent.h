#ifndef STEEPEST_DESCENT_H
#define STEEPEST_DESCENT_H

#include <valarray>
using std::valarray;

valarray<double>  steepest_descent(double (*f)(valarray<double>), valarray<double> (*g)(valarray<double>), valarray<double> x0, double epsilon);

#endif
