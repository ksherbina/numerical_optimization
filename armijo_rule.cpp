/*
Implementation of Inexact Search using Armijo's rule (refer to Chapter 8 in
Bazaraa, M.S., Sherali, H.D., & Shetty, C.M. "Nonlinear Programming (3rd ed)" 
to minimize a function.
*/
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>
#include <tuple>

using std::valarray;

struct Armijo {
    double stepsize;
    int fevals;
};


int comparison (int iter, double (*f)(valarray<double>), const double& derivative, const valarray<double>& xc,
                const valarray<double>& d, double steplength, double theta) {
  double fxc, lc, fxk;
  valarray<double> xk;
  fxc = f(xc);
  lc = fxc + theta * steplength * derivative;
  xk = xc + steplength * d;
  fxk = f(xk);

  int s1 = 18;
  //printf("%d %*.8f %*.8f %*.8f\n", iter, s1, steplength, s1, fxk, s1, lc);
  
  if (fxk <= lc) {
    return 1;
  } else {
    return 0;
  }
}

Armijo armijo_rule(double (*f)(valarray<double>), valarray<double> gradient, valarray<double> x0,
                   valarray<double> direction, double eta, double theta) {
  //User must supply the function f: R^n -> R and the gradient of the function g:R^n.
  Armijo results;
  double alpha, step, descent;
  int counter = 1;
  int function_evaluations = 0;
  
  alpha = 1.0;
  descent = std::inner_product(std::begin(gradient), std::end(gradient), std::begin(direction), 0.0);
  
  int s2 = 10;
  int s3 = 24;

  if (descent < 0) {
    //printf("NOTE: LH is the left hand side of Armijo's rule and RH is the right hand side\n");
    //printf("%s %*s %*s %*s\n","Iteration", s2, "Steplength", s2, "LH", s3, "RH");
    if (comparison(counter, f, descent, x0, direction, alpha, theta) == 1) {
      function_evaluations += 1;
      step = pow(eta, counter) * alpha;
      while (comparison(counter+1, f, descent, x0, direction, step, theta) == 1) {
        function_evaluations += 1;
        alpha = step;
        step = pow(eta, counter) * alpha;
        counter++;
      }
      function_evaluations += 1;
    } else {
      function_evaluations += 1;
      alpha = pow(1.0/eta, counter) * alpha;
      while (comparison(counter+1, f, descent, x0, direction, alpha, theta) == 0) {
        function_evaluations += 1;
        alpha = pow(1.0/eta, counter) * alpha;
        counter++;
      }
      function_evaluations += 1;
    }
    results.stepsize = alpha;
    results.fevals = function_evaluations;
    return results;
    
  } else {
    printf("Current direction is not a descent direction");
    results.stepsize = 0.0;
    results.fevals = 0;
    return results;
  }
}
