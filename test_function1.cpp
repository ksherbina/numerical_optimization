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

double test_func(valarray<double> x) {
    //user-defined function f such that f:R^n->R.
    double sum;
    for (int i = 0; i < (x.size() - 1); i++) {
      sum += ((10.0 - static_cast<double>(i)) * pow((pow(x[i], 2)-x[i+1]), 2));
    }
    
    return pow((x[0]-1.0), 2) + pow((x[9]-1.0), 2) + 10.0 * sum;
}

valarray<double> test_gradient(valarray<double> x) {
    valarray<double> soln(x.size());
    soln[0] = 2.0 * (x[0] - 1.0) + 180.0 * (pow(x[0], 2) - x[1]) * (2.0 * x[0]);
    soln[1] = -180.0 * (pow(x[0], 2) - x[1]) + 160.0 * (pow(x[1], 2) - x[2]) * (2.0 * x[1]);
    soln[2] = -160.0 * (pow(x[1], 2) - x[2]) + 140.0 * (pow(x[2], 2) - x[3]) * (2.0 * x[2]);
    soln[3] = -140.0 * (pow(x[2], 2) - x[3]) + 120.0 * (pow(x[3], 2) - x[4]) * (2.0 * x[3]);
    soln[4] = -120.0 * (pow(x[3], 2) - x[4]) + 100.0 * (pow(x[4], 2) - x[5]) * (2.0 * x[4]);
    soln[5] = -100.0 * (pow(x[4], 2) - x[5]) + 80.0 * (pow(x[5], 2) - x[6]) * (2.0 * x[5]);
    soln[6] = -80.0 * (pow(x[5], 2) - x[6]) + 60.0 * (pow(x[6], 2) - x[7]) * (2.0 * x[6]);
    soln[7] = -60.0 * (pow(x[6], 2) - x[7]) + 40.0 * (pow(x[7], 2) - x[8]) * (2.0 * x[7]);
    soln[8] = -40.0 * (pow(x[7], 2) - x[8]) + 20.0 * (pow(x[8], 2) - x[9]) * (2.0 * x[8]);
    soln[9] = 2.0 * (x[9] - 1.0) - 20.0 * (pow(x[8], 2) - x[9]);
    return soln;
}

valarray<double> test_hessian(valarray<double> x) {
  valarray<double> soln(0.0, x.size() * x.size());
  soln[0] = 2.0 * x[0] + 360.0 * (3.0 * pow(x[0], 2) - x[1]);
  soln[1] = -180.0 * 2.0 * x[0];
  soln[1 * x.size() + 0] = -180.0 * (2.0 * x[0]);  
  soln[1 * x.size() + 1] = -180.0 * x[1] + 320.0 * 3.0 * pow(x[1], 2); 
  soln[1 * x.size() + 2] = -320.0 * x[1]; 
  soln[2 * x.size() + 1] = -160.0 * (2.0 * x[1]);  
  soln[3 * x.size() + 2] = -140.0 * (2.0 * x[2]);
  soln[4 * x.size() + 3] = -120.0 * (2.0 * x[3]);
  soln[5 * x.size() + 4] = -100.0 * (2.0 * x[4]);
  soln[6 * x.size() + 5] = -80.0 * (2.0 * x[5]);
  soln[7 * x.size() + 6] = -60.0 * (2.0 * x[6]);
  soln[8 * x.size() + 7] = -40.0 * (2.0 * x[7]);
  soln[2 * x.size() + 2] = -160.0 * x[2] + 280.0 * 3.0 * pow(x[2], 2);
  soln[3 * x.size() + 3] = -140.0 * x[3] + 240.0 * 3.0 * pow(x[3], 2); 
  soln[4 * x.size() + 4] = -120.0 * x[4] + 200.0 * 3.0 * pow(x[4], 2);
  soln[5 * x.size() + 5] = -100.0 * x[5] + 160.0 * 3.0 * pow(x[5], 2);
  soln[6 * x.size() + 6] = -80.0 * x[6] + 120.0 * 3.0 * pow(x[6], 2);
  soln[7 * x.size() + 7] = -60.0 * x[7] + 80.0 * 3.0 * pow(x[7], 2);
  soln[8 * x.size() + 8] = -40.0 * x[8] + 40.0 * 3.0 * pow(x[8], 2);
  soln[2 * x.size() + 3] = -280.0 * x[2]; 
  soln[3 * x.size() + 4] = -240.0 * x[3]; 
  soln[4 * x.size() + 5] = -200.0 * x[4]; 
  soln[5 * x.size() + 6] = -160.0 * x[5]; 
  soln[6 * x.size() + 7] = -120.0 * x[6]; 
  soln[7 * x.size() + 8] = -80.0 * x[7]; 
  soln[8 * x.size() + 9] = -40.0 * x[8]; 
  soln[9 * x.size() + 8] = -40.0 * x[8];
  soln[9 * x.size() + 0] = 2.0 + 20.0;

  return soln;
}

int main()
{
  valarray<double> sd, modified_newton, testx(0.0, 10);
  testx[0] = -1.2;
 
  for (int i = 0; i < 3; i++) {
    double check;
    check = static_cast<double>(i);
    printf("solution=%2.8f\n", check);
  }
  sd = steepest_descent(test_func, test_gradient, testx, pow(10.0,-3), 200);
  modified_newton = newton(testx, test_func, test_gradient, test_hessian, pow(10.0,-3), "modified", "armijo", 200);

  return 0;
}
