#ifndef NEWTON_1D_H
#define NEWTON_1D_H
 
#include <tuple>

std::tuple<double, double>  newton_1d(double (*f)(double),double (*g)(double),double (*h)(double),double x0, double a, double b, double epsilon, double theta);
 
#endif
