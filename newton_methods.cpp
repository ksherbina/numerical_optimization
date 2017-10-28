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
#include "cholesky.h"
#include "solve_linear_system.h"

using std::valarray;

valarray<double> newton_methods(double (*function)(valarray<double>, int n, double parm),
                      valarray<double> (*gradient)(valarray<double>, int n, double parm),
                      valarray<double> (*hessian)(valarray<double>, int n, double parm),
                      double rosenbrock_parm, valarray<double> x0, int dim, double tolerance,
                      bool modify_diagonal, bool do_linesearch) {
  valarray<double> direction, h, gxk, xk, xc;
  CholeskyFactors hxc;
  double alpha, gradient_norm;
  int counter = 1;
  int max_iters = 50;
  xk = x0;
  gxk = gradient(x0, dim, rosenbrock_parm);
  gradient_norm = std::sqrt(std::inner_product(std::begin(gxk), std::end(gxk), std::begin(gxk), 0.0));
  
  while (counter <= max_iters) {
    xc = xk;
    h = hessian(xc, dim, rosenbrock_parm);
    hxc = cholesky(h, dim, modify_diagonal);
    if (modify_diagonal) {
      direction = solve_linear_system(hxc,-1*gxk,dim);
    } else {
      // Must check if Hessian is positive definite. If not, take
      // steepest descent direction rather than newton direction
      std::cout<<"Hessian may not be positive definite"<<std::endl;
    }
    if (do_linesearch) {
      alpha = 2.0;
      //alpha = armijo_rule(function, gradient, xc, direction, 2.0, 0.5);
    } else {
      alpha = 1.0;
    }
    xk = xc + alpha * direction;
    gxk = gradient(xk, dim, rosenbrock_parm);
    gradient_norm = std::sqrt(std::inner_product(std::begin(gxk), std::end(gxk), std::begin(gxk), 0.0));
    counter++;
    if (counter > max_iters) {
      std::cout<<"Exceeded maximum number of iterations!"<<std::endl;
      break;
    }

  }
  return xk;
}
