/*
Implementation of Inexact Search using Armijo's rule (refer to Chapter 8 in
Bazaraa, M.S., Sherali, H.D., & Shetty, C.M. "Nonlinear Programming (3rd ed)" 
to minimize a function.
*/
#include <iostream>
#include <cmath>
#include <stdio.h> //include printfin
#include <valarray>
#include <typeinfo>

using std::valarray;

double armijo_rule(int n,double (*f)(valarray<double>),valarray<double> (*g)(valarray<double>),valarray<double> xk, valarray<double> d, double a, double eta, double theta)
{
  //User must supply the function f: R^n -> R and the gradient of the function g:R^n.
  double alpha;
  alpha=a;
 
  return alpha;
}
