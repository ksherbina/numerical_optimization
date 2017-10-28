/*
Implementation of either Pure or Modified Newton method to find a 
direction.
*/
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>
#include "cholesky.h"
#include "solve_linear_system.h"

using std::valarray;

valarray<double> newton_direction(double (*function)(valarray<double>),valarray<double> (*gradient)(valarray<double>),valarray<double> (*hessian)(valarray<double>),valarray<double> x0, int dim, double epsilon, bool modify_diagonal, bool do_linesearch)
{
  valarray<double> d, h, gxk, xk, xc;
  CholeskyFactors hxc;
  double alpha;
  valarray<double> eps (epsilon,dim);   
  int counter=1;
  xk=x0;
  gxk=gradient(xk);
  
  while (std::abs(gxk)>eps) {
    xc=xk;
    h=hessian(xc);
    hxc=cholesky(h, dim, modify_diagonal);
    if (modify_diagonal) {
      d = solve_linear_system(hxc,-1*gxk,dim);
    } else {
    // Must check if Hessian is positive definite. If not, take 
    // steepest descent direction rather than newton direction
    }
    if (linesearch) {
      alpha=armijo_rule(function,gradient,xc,d,2.0,0.5);
    } else {
      alpha=1.0;
    }
    xk=xc+alpha*d;
    gxk=gradient(xk);
    counter+=1;
  }
  return d;
}
