#ifndef ARMIJO_RULE_H
#define ARMIJO_RULE_H
 
#include <valarray>
using std::valarray;

double armijo_rule(int n,double (*f)(valarray<double>),valarray<double> (*g)(valarray<double>),valarray<double> xk, valarray<double> d, double a, double eta, double theta);
 
#endif
