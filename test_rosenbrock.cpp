//Function, inputs, and parameters for problem 1 of homework 2a in ISE 520 Fall 2017
#include <limits>
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "cholesky.h"
#include "steepest_descent.h"
#include "golden_section_search.h"
#include "solve_linear_system.h"
#include "newton.h"
#include "armijo_rule.h"

using std::valarray;
using std::string;

double rosenbrock(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  double value = 0.0;
  int n = x.size();
  double alpha = 100.0;
  double first_term, second_term;
  for (int i = 0; i < n / 2; i++) {
    first_term = (x[2 * i + 1]-(x[2 * i] * x[2 * i]))*(x[2 * i + 1]-(x[2 * i] * x[2 * i]));
    second_term = (1.0 - x[2 * i]) * (1.0 - x[2 * i]);
    value = value + (alpha * first_term + second_term);
  }
  return value;
}

valarray<double> rosenbrock_gradient(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  double alpha = 100;
  valarray<double> value(0.0, n);
  for (int i = 0; i < n / 2; i++) {
    value[0] = -4.0 * alpha * (x[2 * i + 1] - (x[2 * i] * x[2 * i])) * x[2 * i] - 2.0 * (1 - x[2 * i]) + value[0];
    value[1] = x[2 * i + 1] - (x[2 * i] * x[2 * i]) + value[1];
    std::cout<<value[1]<<std::endl;
  }
  value[1] = 2.0 * alpha * value[1];
  return value;
}

valarray<double> rosenbrock_hessian(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  double alpha = 100;
  valarray<double> value(0.0, n);
  for (int i = 0; i < n / 2; i++) {
    value[0] = 4.0 * alpha * (x[2 * i + 1] - 3.0 * (x[2 * i] * x[2 * i])) + 2.0 + value[0];
    value[1] = 4.0 * alpha * x[2 * i + 1] + value[1];
    value[2] = 4.0 * alpha * x[2 * i] + value[2];
  }
  value[3] = -2.0 * alpha;
  return value;
}

int main()
{
  valarray<double> sd, modified_newton;
  valarray<double> testx = { -1.2, -1.2 };
 
  sd = steepest_descent(rosenbrock, rosenbrock_gradient, testx, pow(10.0,-2), 200);
  modified_newton = newton(testx, rosenbrock, rosenbrock_gradient, rosenbrock_hessian, pow(10.0,-2), "modified", "armijo", 200);
  return 0;
}
