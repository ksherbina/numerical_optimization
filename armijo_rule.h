#ifndef ARMIJO_RULE_H
#define ARMIJO_RULE_H
 
#include <valarray>
using std::valarray;

struct Armijo {
    double stepsize;
    int fevals;
};

Armijo armijo_rule(double (*f)(valarray<double>), valarray<double> gradient, valarray<double> x0,
                   valarray<double> direction, double eta, double theta);
 
#endif
