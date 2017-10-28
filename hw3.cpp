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
#include "newton_methods.h"

double rosenbrock(valarray<double> x, int n, double alpha) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  double value = 0.0;
  double first_term, second_term;
  for (int i = 0; i < n / 2; i++) {
    first_term = (x[2 * i + 1]-(x[2 * i] * x[2 * i]))*(x[2 * i + 1]-(x[2 * i] * x[2 * i]));
    second_term = (1.0 - x[2 * i]) * (1.0 - x[2 * i]);
    value = value + (alpha * first_term + second_term);
  }
  return value;
}

valarray<double> rosenbrock_gradient(valarray<double> x, int n, double alpha) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  valarray<double> value(0.0, n);
  for (int i = 0; i < n / 2; i++) {
    value[0] = -4.0 * alpha * (x[2 * i + 1] - (x[2 * i] * x[2 * i])) * x[2 * i] - 2.0 * (1 - x[2 * i]) + value[0];
    value[1] = x[2 * i + 1] - (x[2 * i] * x[2 * i]) + value[1];
    std::cout<<value[1]<<std::endl;
  }
  value[1] = 2.0 * alpha * value[1];
  return value;
}

valarray<double> rosenbrock_hessian(valarray<double> x, int n, double alpha) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  valarray<double> value(0.0, n*n);
  for (int i = 0; i < n / 2; i++) {
    value[0] = 4.0 * alpha * (x[2 * i + 1] - 3.0 * (x[2 * i] * x[2 * i])) + 2.0 + value[0];
    value[1] = 4.0 * alpha * x[2 * i + 1] + value[1];
    value[2] = 4.0 * alpha * x[2 * i] + value[2];
  }
  value[3] = -2.0 * alpha;
  return value;
}


int main() {
  int ncol = 3;
  CholeskyFactors check_chol;
  double check_function;
  valarray<double> check_solv, check_derivative, check_hessian, check_newton;
  valarray<double> M = { 2.0, -1.0, 1.0, -1.0, 3.0, 0.0, 1.0, 0.0, 5.0 };
  valarray<double> b = { 1.0, -2.0, 3.0 };
  valarray<double> initial_point = { -1.0, -1.0 };
  double epsilon = pow(10.0, -6);
  check_chol = cholesky(M, ncol, 0);
  check_solv = solve_linear_system(check_chol, b, ncol);
  check_function = rosenbrock(initial_point, 2, 100.0);
  std::cout<<"Initial Rosenbrock function value = "<<check_function<<std::endl;
  check_derivative = rosenbrock_gradient(initial_point, 2, 100.0);
  for (int i = 0; i < check_derivative.size(); i++) {
    printf("f'[%d] = %.8f \n", i, check_derivative[i]);
  }
  check_hessian = rosenbrock_hessian(initial_point, 2, 100.0);
  for (int i = 0; i < check_hessian.size(); i++) {
    printf("f''[%d] = %.8f \n", i, check_hessian[i]);
  }
  
  check_newton = newton_methods(rosenbrock, rosenbrock_gradient, rosenbrock_hessian,
                                initial_point, ncol, epsilon, 0, 0);
  
  std::cout<<std::numeric_limits<double>::epsilon()<<std::endl;
  
  return 0;
}

