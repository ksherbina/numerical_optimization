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
#include "damped_bfgs.h"
#include "solve_qp_froth.h"
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
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = 0.0;
  value[1] = 0.0;
  value[2] = 0.0;
  value[3] = 12.0;
  return value;
}

valarray<double> constraints(valarray<double> x) {
  int n = x.size();
  valarray<double> value(0.0, n);
  value[0] = x[0] - x[1] * (2.0 - x[1] * (5.0 - x[1])) - 13.0 - 50.0;
  value[1] = x[0] - x[1] * (14.0 - x[1] * (1.0 + x[1])) - 29.0 - 50.0;
  return value;
}

valarray<double> c1_gradient(valarray<double> x) {
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
  valarray<double> modified_newton, bfgs, L, x;
  valarray<double> testx = {0.5, -2.0};
  Armijo steps;

  modified_newton = newton(testx, froth, froth_gradient, froth_hessian, pow(10.0,-6), "modified", "armijo", 200);
  bfgs = newton(testx, froth, froth_gradient, froth_hessian, pow(10.0,-6), "bfgs", "armijo", 200);

  QpSolutions frothqp;
  valarray<double> L0(0.0, testx.size() * testx.size()), plambda;
  valarray<double> lambda(1.0, testx.size()); // Initialize Lagrange multipliers
  double gradient_norm;

  // Evaluate function, gradient, and constraints at current point
  fx = froth(testx);
  gx = froth_gradient(testx);
  hx = froth_hessian(textx);
  cx = constraints(testx);
  cx1g = c1_gradient(testx);
  cx2g = c2_gradient(testx);
  c1xh = c1_hessian(testx);
  c2xh = c2_hessian(testx);

  gradient_norm = std::sqrt(std::inner_product(std::begin(gx), std::end(gx), std::begin(gx), 0.0));

  for (int i = 0; i < L0.size(); i += (x.size() + 1)) {
    L0[i] = 1.0;
  } 
 
  while (gradient_norm > epsilon) { 
     frothqp = solve_qp_froth(fx, gx, hx, L0, cx, cx1g, cx2g, c1xh, c2xh);
     plambda = lambda - frothqp.multipliers;
     // Find mu_k as in equation 18.36 in Nocedal & Wright with sigma = 1;
     // Do Armijo's rule to find stepsize alpha where gx0 is the directional derivative of the function phi as 
     // described in Algorithm 18.3 in Nocedal & Wright
     steps = armijo_rule(func, gx0, x0, frothqp.direction, 1.0, 0.5);
     func_evals += steps.fevals;
     x = testx + steps.stepsize * frothqp.direction;
     delta = x - textx;
     gamma = L - L0;
     fx = f(x);
     gx = g(x);
     // Reevaluate all constraints and their derivatives
     L = damped_bfgs(L0, delta, gamma)
     testx = x;
     L0 = L;
  }
  
  return 0;
}
