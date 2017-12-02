//Implementation of steepest descent to find best direction in which to perform 
//minimize a function.

#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>
#include <tuple>
#include "armijo_rule.h"
#include "golden_section_search.h"

using std::valarray;

valarray<double> steepest_descent(double (*f)(valarray<double>), valarray<double> (*g)(valarray<double>), valarray<double> x0, double epsilon) {

  int k = 0;  
  double fx, stepsize, gnorm, a, b;
  valarray<double> x, gx, d, newx;
  
  x = x0;
  fx = f(x);
  gx = g(x);
  gnorm = std::sqrt(std::inner_product(std::begin(gx), std::end(gx), std::begin(gx), 0.0));
  
  std::cout << "Enter left endpoint for Golden Section Search: ";
  std::cin >> a;
  std::cout << "Enter right endpoint for Golden Section Search: ";
  std::cin >> b;

  int s1=22,s2=16,s3=14;
  printf("%s %s %*s %*s %*s \n","Iteration (i)", "x_i[0]", s2, "x_i[1]", s2, "f(x_i)", s2, "norm of gradient");

  while (gnorm > epsilon) {
    printf("%d %*.8f %*.8f %*.8f %*.8f \n", k, s1, x[0], s3, x[1], s3, fx, s3, gnorm);
    d = -gx;
    newx = golden_section_search(x.size(), f, x, a, b, d, pow(10.0, -3));
    x = newx;
    fx = f(x);
    gx = g(x);
    gnorm = std::sqrt(std::inner_product(std::begin(gx), std::end(gx), std::begin(gx), 0.0));
    k++;
    if (k > 20) {
      printf("Maximum iterations exceeded");
      break;
    }
  }

  return x;
}

