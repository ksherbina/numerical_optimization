#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include "cholesky.h"
#include "solve_linear_system.h"
#include "armijo_rule.h"

using std::valarray;
using std::string;

double scalar_product_of_two_vectors(valarray<double> u, valarray<double> v) {
  /*
  Multiply a 1xn vector u by a nx1 vector v 
  */
  double sum = 0.0;
  for (int j = 0; j < u.size(); j++) {
    sum += (u[j] * v[j]);
  }
  return sum;
}

valarray<double> newton(valarray<double> x0, double (*func)(valarray<double>), valarray<double> (*grad)(valarray<double>),
                                valarray<double> (*hess)(valarray<double>), double epsilon, string newton_type, string linsearch_type, int MAX_ITER) {

  valarray<double> x, gx, gx0, hx, d;
  CholeskyFactors hxc;
  Armijo steps;
  double gradient_norm, stepsize, fx;

  int func_evals = 0;
  int counter = 0;

  x = x0;
  gx = grad(x);

  int gradient_evaluations = 1;
  
  gradient_norm = std::sqrt(std::inner_product(std::begin(gx), std::end(gx), std::begin(gx), 0.0)); 
  /*
  hx = hess(x);
  hxc = cholesky(hx, x.size());
  d = solve_linear_system(hxc, -gx, gx.size());
  std::cout<<"Check direction: "<<scalar_product_of_two_vectors(gx, d)<<std::endl;
  */

  int s1=22,s2=16,s3=14;
  printf("%s %s %*s %*s %*s \n","Iteration (i)", "x_i[0]", s2, "x_i[1]", s2, "f(x_i)", s2, "norm of gradient");

  while (gradient_norm > epsilon) {
    fx = func(x);
    printf("%d %*.8f %*.8f %*.8f %*.8f \n", counter, s1, x[0], s3, x[1], s3, fx, s3, gradient_norm);
    hx = hess(x);
    hxc = cholesky(hx, x.size());
    d = solve_linear_system(hxc, -gx, gx.size());
    /*
    for (int i = 0; i < gx.size(); i++) {
      printf("gx=%2.8f\n", gx[i]);
    }
    for (int i = 0; i < d.size(); i++) {
      printf("d=%2.8f\n", d[i]);
    }
    */
    if (scalar_product_of_two_vectors(gx, d) < 0.0) {
      x0 = x;
      gx0 = gx;
      if (newton_type == "pure") {
        x = x0 + d;
      } else {
        steps = armijo_rule(func, gx0, x0, d, 2.0, 0.5);
        func_evals += steps.fevals;
        x = x0 + steps.stepsize * d;
      }
      gx = grad(x);
      gradient_evaluations += 1;
      gradient_norm = std::sqrt(std::inner_product(std::begin(gx), std::end(gx), std::begin(gx), 0.0));
    } else {
      std::cout<<"Direction is not a direction of descent"<<std::endl;
      break;
    }
    counter += 1;
    if (counter > MAX_ITER) {
      printf("Maximum iterations exceeded \n");
      break;
    }
  }
  std::cout<<func_evals<<" function evaluations"<<std::endl;
  std::cout<<gradient_evaluations<<" gradient evaluations"<<std::endl;
  return x;
}
