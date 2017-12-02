//Function, inputs, and parameters for problem 1 of homework 2a in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "steepest_descent.h"
#include "golden_section_search.h"

using std::valarray;

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


int main()
{
  /*
  int n=3;
  valarray<double> solv1;
  valarray<double> M {2.0,-1.0,1.0,-1.0,3.0,0.0,1.0,0.0,5.0};
  CholeskyFactors chol1;
  chol1=cholesky(M,n);
  valarray<double> b {1.0,-2.0,3.0};
  solv1 = solve_linear_system(chol1,b,n);
  for (int i=0;i<n;i++) {
      printf("solution=%2.8f\n",solv1[i]);
  }
  */

  valarray<double> testx, check_sd;
  testx = {-4.0, -2.0};
  double lx = -1.0;
  double rx = 2.0;
  check_sd = steepest_descent(testf1, test_gradient1, testx, pow(10.0,-3), 50);
  
  /*
  double a = 0.4;
  double b = 1.5;
  double d = 
  gss=golden_section_search(n,func1,x,a,b,d,error);
  */
  return 0;
}
