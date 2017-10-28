#ifndef NEWTON_METHODS_H
#define NEWTON_METHODS_H
 
#include <valarray>
using std::valarray;

double newton_methods(double (*function)(valarray<double>, int n, double parm),
                      valarray<double> (*gradient)(valarray<double>, int n, double parm),
                      valarray<double> (*hessian)(valarray<double>, int n, double parm),
                      valarray<double> x0, int dim, double tolerance, bool modify_diagonal,
                      bool do_linesearch);
 
#endif