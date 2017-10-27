//Function, inputs, and parameters for problem 1 of homework 2a in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "armijo_rule.h"

using std::valarray;

double func1b(valarray<double> x) {
  //user-defined function f such that f:R^n->R.
  return 2*x[0]*x[1]+(10/x[1])+(5*x[0]);
}

valarray<double> grad1b(valarray<double> x) {
  //gradient of func1a; in R^n
  valarray<double> g (2);
  g[0]=2*x[1]-(5/pow(x[0],2));
  g[1]=2*x[0]-(10/pow(x[1],2));
  return g;
}

int main()
{
  clock_t t;
  int n=2; //User must specify number of dimensions
  
  valarray<double> d(n), x0(n), ar1(n), ar2(n), ar3(n);
  double stepsize=1;
  double eta=2.0;
  double theta=0.5;
  
  //TA's Problem 1b
  d={1.0,1.0};
  x0={1.0,1.0};
  std::cout<<"Minimize objective function in problem 1b"<<std::endl;
  t=clock();
  std::cout<<"Running Armijio's Rule Inexact Line Search..."<<std::endl;
  ar2=armijo_rule(func1b,grad1b,x0,d,stepsize,eta,theta);
  printf("Final output = [%2.8f %2.8f]\n",ar2[0],ar2[1]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  return 0;
}
