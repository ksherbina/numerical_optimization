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

using std::valarray;

double armijo_rule(int n,double (*f)(valarray<double>),valarray<double> gxk,valarray<double> xk, valarray<double> d, double a, double eta, double theta)
{
  //User must supply the function f: R^n -> R and the gradient of the function g:R^n.
  valarray<double> xc(n);
  double alpha,fxc,lc,prod;
  
  int t=1;
  alpha=a;
  double fxk=f(xk);
  xc=xk+alpha*d;
  fxc=f(xc);
  prod=std::inner_product(std::begin(gxk),std::end(gxk),std::begin(d),0.0);
  std::cout<<prod<<std::endl;
 
  return alpha;
}
