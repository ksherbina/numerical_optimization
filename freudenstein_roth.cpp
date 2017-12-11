// Freudenstein and Roth function
#include <limits>
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "cholesky.h"
#include "solve_linear_system.h"
#include "newton.h"
#include "quasi_newton.h"
#include "armijo_rule.h"

using std::valarray;
using std::string;

double froth(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  double first_term, second_term;
  first_term = x[0] - x[1] * (2.0 - x[1] * (5.0 - x[1])) - 13.0;
  second_term = x[0] - x[1] * (14.0 - x[1] * (1.0 + x[1])) - 29.0;
  return first_term + second_term;
}

valarray<double> froth_gradient(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 2.0;
  value[1] = 12.0 * x[1] - 16.0;
  return value;
}

valarray<double> froth_hessian(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  value[3] = 12.0;
  return value;
}

valarray<double> constraints(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = x[0] - x[1] * (2.0 - x[1] * (5.0 - x[1])) - 13.0 - 30.0;
  value[1] = x[0] - x[1] * (14.0 - x[1] * (1.0 + x[1])) - 29.0 - 30.0;
  return value;
}

valarray<double> c1_gradient(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 1.0;
  value[1] = -3.0 * pow(x[1], 2) + 10.0 * x[1] - 2.0;
  return value;
}

valarray<double> c2_gradient(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 1.0;
  value[1] = 3.0 * pow(x[1], 2) + 2.0 * x[1] - 14.0;
  return value;
}

valarray<double> c1_hessian(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  value[3] = -6.0 * x[1] + 10.0;
  return value;
}

valarray<double> c2_hessian(valarray<double> x) {
  //n must be even
  //alpha is a user-specified parameter (ex. 1 or 100)
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  value[3] = 6.0 * x[1] + 2.0;
  return value;
}

int main()
{
  valarray<double> modified_newton, bfgs;
  valarray<double> testx = {0.5, -2.0};
 
  //modified_newton = newton(testx, rosenbrock, rosenbrock_gradient, rosenbrock_hessian, pow(10.0,-2), "modified", "armijo", 200);
  return 0;
}
