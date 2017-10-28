/*
 Solutions to Homework 3 from ISE 520 Fall 2017 by Katrina Sherbina.
 */
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>
#include "cholesky.h"
#include "solve_linear_system.h"

double rosenbrock(valarray<double> x, int n, double alpha) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  double value = 0.0;
  double first_term;
  double second_term;
  for (int i = 0; i < n / 2; i++) {
    first_term = (x[2 * i + 1]-(x[2 * i] * x[2 * i]))*(x[2 * i + 1]-(x[2 * i] * x[2 * i]));
    second_term = (1.0 - x[2 * i]) * (1.0 - x[2 * i]);
    value = value + (alpha * first_term + second_term);
  }
  return value;
}

valarray<double> grad1a(valarray<double> x) {
  //gradient of func1a; in R^n
  valarray<double> g (2);
  g[0]=8.0*x[0]+4.0*x[1];
  g[1]=4.0*x[0]+8.0*x[1]-12.0;
  return g;
}


int main() {
  int ncol = 3;
  CholeskyFactors check_chol;
  double check_function;
  valarray<double> check_solv, check_newton;
  valarray<double> M = { 2.0, -1.0, 1.0, -1.0, 3.0, 0.0, 1.0, 0.0, 5.0 };
  valarray<double> b = { 1.0, -2.0, 3.0 };
  valarray<double> initial_point = { -1.0, -1.0 };
  
  check_chol = cholesky(M, ncol, 0);
  check_solv = solve_linear_system(check_chol, b, ncol);
  check_function = rosenbrock(initial_point, 2, 100.0);
  std::cout<<"Initial Rosenbrock function value = "<<check_function<<std::endl;
  /*
  check_newton = newton_methods(double (*function)(valarray<double>),valarray<double> (*gradient)(valarray<double>),
                                valarray<double> (*hessian)(valarray<double>), valarray<double> x0, int dim,
                                double tolerance, bool modify_diagonal, bool do_linesearch);
  */
  std::cout<<std::numeric_limits<double>::epsilon()<<std::endl;
  
  return 0;
}

