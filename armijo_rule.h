#ifndef ARMIJO_RULE_H
#define ARMIJO_RULE_H
 
#include <valarray>
using std::valarray;

double armijo_rule(double (*f)(valarray<double>),valarray<double> (*g)(valarray<double>),valarray<double> x0,valarray<double> d, double eta, double theta);
 
#endif
