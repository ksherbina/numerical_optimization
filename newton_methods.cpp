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

valarray<double> newton_methods(double tolerance) {
  valarray<double> direction, h, gxk, xk(0.0,2), xc;
  double alpha, gradient_norm;
  int counter = 1;
  int max_iters = 50;
  
  return xk;
}
