/*
Find Cholesky factors of a square matrix and check if the 
matrix is positive definite.
*/
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm> //for std::max_element

using std::valarray;

valarray<double> newton_methods(double (*function)(valarray<double>, int n, double parm),
                      valarray<double> (*gradient)(valarray<double>, int n, double parm),
                      valarray<double> (*hessian)(valarray<double>, int n, double parm),
                      valarray<double> x0, int dim, double tolerance, bool modify_diagonal,
                      bool do_linesearch) {
  valarray<double> direction, h, gxk, xc;
  double alpha, gradient_norm;
  int counter = 1;
  int max_iters = 50;
  
  valarray<double> xk(0.0,2);
  
  return xk;
}
