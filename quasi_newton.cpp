// BFGS and DFP update of Hessian matrix
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include "cholesky.h"
#include "newton.h"
#include "cholesky_ranktwo.h"
#include "solve_linear_system.h"
#include "armijo_rule.h"

using std::valarray;
using std::string;

double scalar_product_of_2_vectors(valarray<double> u, valarray<double> v) {
  /*
  Multiply a 1xn vector u by a nx1 vector v 
  */
  double sum = 0.0;
  for (int j = 0; j < u.size(); j++) {
    sum += (u[j] * v[j]);
  }
  return sum;
}

valarray<double> multiply_matrix_by_vector(valarray<double> M, valarray<double> z) {
  /*
  Multiply a nxn matrix M by a nx1 vector z 
  */
  valarray<double> prod(0.0, z.size());
  double sum;
  for (int j = 0; j < M.size(); j += z.size()) {
    sum = 0.0;
    for (int k = 0; k < z.size(); k++) {
      sum += (M[j + k] * z[k]);
    }
    prod[j/z.size()] = sum;
  }
  return prod;
}

valarray<double> multiply_vector_by_its_transpose(valarray<double> y) {
  valarray<double> mat(0.0, y.size() * y.size());
  for (int j = 0; j < y.size(); j++) {
    for (int k = 0; k < y.size(); k++) {
      mat[j*y.size()+k] = y[j] * y[k];
    }
  }
  return mat;
}

valarray<double> quasi_newton(valarray<double> x0, double (*func)(valarray<double>), valarray<double> (*grad)(valarray<double>),
                                double epsilon, string newton_type, string linsearch_type, int MAX_ITER) {

  valarray<double> x, gx, gx0, delta, gamma, d, H0(0.0, x0.size() * x0.size()), H, term1, term2, hgamma, hdelta, add_rank1, subtract_rank1;
  CholeskyFactors hxc0, hxc;
  Armijo steps;
  double fx, gradient_norm, curvature, sum, stepsize;
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
  
  if (newton_type == "dfp") {
    for (int i = 0; i < H0.size(); i += (x.size() + 1)) {
      H0[i] = 1.0;
    } 
    d = multiply_matrix_by_vector(-H, gx);
  }

  if (newton_type == "bfgs") {
    for (int i = 0; i < H0.size(); i += (x.size() + 1)) {
      H0[i] = 1.0;
    } 
    hxc = cholesky(H0, x.size());
    d = solve_linear_system(hxc, -gx, gx.size());
  }

  int s1=22,s2=16,s3=14;
  printf("%s %s %*s %*s %*s \n","Iteration (i)", "x_i[0]", s2, "x_i[1]", s2, "f(x_i)", s2, "norm of gradient");

  while (gradient_norm > epsilon) {
    fx = func(x);
    printf("%d %*.8f %*.8f %*.8f %*.8f \n", counter, s1, x[0], s3, x[1], s3, fx, s3, gradient_norm);
    if (counter > 0) {
      if (newton_type == "dfp") {
        delta = x - x0;
        gamma = gx - gx0;
        hgamma = multiply_matrix_by_vector(H0, gamma);
        term1 = multiply_vector_by_its_transpose(hgamma)/scalar_product_of_2_vectors(gamma, hgamma);
        term2 = multiply_vector_by_its_transpose(delta)/scalar_product_of_2_vectors(gamma, delta);
        H = H0 - term1 - term2;
        d = multiply_matrix_by_vector(-H, gx);
        H0 = H;
      }
      if (newton_type == "bfgs") {   
        delta = x - x0;
        gamma = gx - gx0;
        curvature = std::inner_product(std::begin(delta), std::end(delta), std::begin(gamma), 0.0);
        add_rank1 = multiply_vector_by_its_transpose(gamma) / scalar_product_of_2_vectors(gamma, delta);
        hdelta = multiply_matrix_by_vector(H0, delta);
        subtract_rank1 = hdelta / scalar_product_of_2_vectors(delta, hdelta);
        if (curvature > 0) {
          hxc = cholesky_ranktwo(add_rank1, subtract_rank1, hxc0);
          d = solve_linear_system(hxc, -gx, gx.size());
          H0 = H;
          hxc0 = hxc;
        } 
      } 
    } 
    if (scalar_product_of_2_vectors(gx, -d) < 0) {
      x0 = x;
      gx0 = gx;
      steps = armijo_rule(func, gx0, x0, d, 2.0, 0.5);
      func_evals += steps.fevals;
      x = x0 + steps.stepsize * d;
      gx = grad(x);
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
