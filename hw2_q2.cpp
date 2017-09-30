//Function, inputs, and parameters for problem 1 of homework 2a in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include "newton_1d.h"

double func2a(double x) {
    //user-defined function f such that f:R->R.
    return pow(x,3)-3*pow(x,2);
}

double grad2a(double x) {
  //gradient of func1a; in R^n
  return 3.0*pow(x,2)-6.0*x;
}

double hessian2a(double x) {
  //gradient of func1a; in R^n
  return 6.0*x-6.0;
}

int main()
{
  clock_t t;
  double x0,a,b,xn,fn;
  
  double epsilon=pow(10.0,-8);
  double theta=0.5;

  //Problem 1a
  x0=1.2;
  a=1.0;
  b=4.0;

  /*
  double f1=func2a(x0);
  printf("f1=%2.8f\n",f1);
  double g1=grad2a(x0);
  printf("g1=%2.8f\n",g1);
  double h1=hessian2a(x0);
  printf("h1=%2.8f\n",h1);
   */

  std::cout<<"Minimize f(x,y)=x^3-3*x^2"<<std::endl;
  t=clock();
  std::cout<<"Running 1-Dimensional Newton's Method..."<<std::endl;
  
  std::tie(xn, fn)=newton_1d(func2a,grad2a,hessian2a,x0,a,b,epsilon,theta);
  printf("The minimum of the function over [%2.8f, %2.8f] is %2.8f and occurs at %2.8f\n",a,b,fn,xn);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  
  return 0;
}
