#ifndef NEWTON_LINESEARCH_H
#define NEWTON_LINESEARCH_H

#include <tuple>

std::tuple<double, double>  newton_linesearch(double (*f)(double),double (*g)(double),double (*h)(double),double x0, double a, double b, double epsilon, double theta);

#endif

