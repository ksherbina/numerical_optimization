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

int main() {
  int ncol = 3;
  CholeskyFactors check_chol;
  valarray<double> check_solv;
  valarray<double> M = { 2.0, -1.0, 1.0, -1.0, 3.0, 0.0, 1.0, 0.0, 5.0 };
  valarray<double> b = { 1.0, -2.0, 3.0 };
  
  check_chol = cholesky(M, ncol, 0);
  check_solv = solve_linear_system(check_chol, b, ncol);
  std::cout<<std::numeric_limits<double>::epsilon()<<std::endl;
  
  return 0;
}

