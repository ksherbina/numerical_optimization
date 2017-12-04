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
#include "quasi_newton.h"
#include "armijo_rule.h"

using std::valarray;
using std::string;

double testf1(valarray<double> x) {
    //user-defined function f such that f:R^n->R.
    return pow(x[0], 2) + x[0] * x[1] - 3.0 * sin(x[1]) + (pow(x[1], 2) / 5.0);
}

valarray<double> test_gradient1(valarray<double> x) {
    valarray<double> soln(x.size());
    soln[0] = 2.0 * x[0] + x[1];
    soln[1] = x[0] + (2.0 * x[1] / 5.0) - 3.0 * cos(x[1]);
    return soln;
}

valarray<double> test_hessian1(valarray<double> x) {
  valarray<double> soln(x.size() * x.size());
  soln[0] = x[1];
  soln[1] = 2.0 * x[0];
  soln[2] = 1.0;
  soln[3] = (2.0 / 5.0) + 3.0*sin(x[1]);
}


int main()
{
 
  int n=3;
  double sum, modified_newton_eval, sd_eval;
  valarray<double> solv1, check, testx, modified_newton, dfp, bfgs;
  valarray<double> M {2.0,-1.0,1.0,-1.0,3.0,0.0,1.0,0.0,5.0};
  CholeskyFactors chol1, chol2;
  chol1=cholesky(M,n);
  valarray<double> b {1.0,-2.0,3.0};
  solv1 = solve_linear_system(chol1,b,n);
  /*
  std::cout<<"Check solve_linear_system()"<<std::endl;
  for (int i=0;i<n;i++) {
      printf("solution=%2.8f\n",solv1[i]);
  }
  
  std::cout<<"Cholesky factors of identity matrix:"<<std::endl;
  valarray<double> I {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  chol2=cholesky(I,n);

  std::cout<<"Check matrix by vector mulitplication:"<<std::endl;
  valarray<double> prod(0.0, solv1.size());
  for (int j = 0; j < M.size(); j += solv1.size()) {
    sum = 0.0;
    for (int k = 0; k < solv1.size(); k++) {
      std::cout<<"j = "<<j<<" and k = "<<k<<std::endl;
      sum += (M[j + k] * solv1[k]);
    }
    prod[j/solv1.size()] = sum;
  }
  for (int i = 0; i < prod.size(); i++) {
    printf("solution=%2.8f\n", prod[i]);
  }

  std::cout<<"Check multiply_vector_by_its_transpose:"<<std::endl;
  double val;
  for (int j = 0; j < b.size(); j++) {
    for (int k = 0; k < b.size(); k++) {
      std::cout<<"j = "<<j<<" and k = "<<k<<std::endl;
      val = b[j] * b[k];
      printf("solution=%2.8f\n", val);
    }
  }
  std::cout<<"Check scalar_product_of_two_vectors:"<<std::endl;
  double sum2 = 0.0;
  for (int j = 0; j < b.size(); j++) {
    sum2 += (b[j] * b[j]);
    printf("intermediate=%2.8f\n", sum2);
  }
  printf("solution=%2.8f\n", sum2);

  std::cout<<"Check valarray subtraction:"<<std::endl;
  check = b - b;
  for (int j = 0; j < b.size(); j++) {
    printf("b-b=%2.8f\n", check[j]);
  }
  */
  std::cout<<std::numeric_limits<double>::epsilon()<<std::endl;

  std::cout<<"Check modified newton:"<<std::endl;
  testx = {-4.0, -2.0};
  modified_newton = newton(testx, testf1, test_gradient1, test_hessian1, pow(10.0,-3), "modified", "armijo", 50);
  modified_newton_eval = testf1(modified_newton);
  std::cout<<"Minimum = "<<modified_newton_eval<<std::endl;
  std::cout<<"which occurs at the point "<<std::endl;
  for (int j = 0; j < modified_newton.size(); j++) {
    std::cout<<modified_newton[j]<<std::endl;
  }

  std::cout<<"Check DFP:"<<std::endl;
  dfp = quasi_newton(testx, testf1, test_gradient1, pow(10.0,-3), "dfp", "armijo", 50);
  //std::cout<<"Check BFGS:"<<std::endl;
  //bfgs = quasi_newton(testx, testf1, test_gradient1, pow(10.0,-3), "bfgs", "armijo", 50);
  /*
  dfp_eval = testf1(dfp_newton);
  std::cout<<"Minimum = "<<dfp_eval<<std::endl;
  std::cout<<"which occurs at the point "<<std::endl;
  for (int j = 0; j < dfp.size(); j++) {
    std::cout<<dfp[j]<<std::endl;
  }
  */
  /*
  valarray<double> check_sd;
  double lx = -1.0;
  double rx = 2.0;
  check_sd = steepest_descent(testf1, test_gradient1, testx, pow(10.0,-3), 50);
  sd_eval = testf1(check_sd);
  std::cout<<"Minimum = "<<sd_eval<<std::endl;
  std::cout<<"which occurs at the point "<<std::endl;
  for (int j = 0; j < check_sd.size(); j++) {
    std::cout<<check_sd[j]<<std::endl;
  }
  */
  return 0;
}
