//Function, inputs, and parameters for problem 1 of homework 2a in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "armijo_rule.h"

using std::valarray;

double func1a(valarray<double> x) {
    //user-defined function f such that f:R^n->R.
    return -12.0*x[1]+4.0*(pow(x[0],2.0))+4.0*(pow(x[1],2.0))+4*x[0]*x[1];
}

valarray<double> grad1a(valarray<double> x) {
    //gradient of func1a; in R^n
    valarray<double> g (2);
    g[0]=8.0*x[0]+4.0*x[1];
    g[1]=4.0*x[0]+8.0*x[1]-12.0;
    return g;
}

double func1b(valarray<double> x) {
  //user-defined function f such that f:R^n->R.
  return 2*x[0]*x[1]+(10/x[0])+(5*x[1]/pow(x[0],2));
}

valarray<double> grad1b(valarray<double> x) {
  //gradient of func1a; in R^n
  valarray<double> g (2);
  g[0]=2*x[1]-(10/pow(x[0],2))-(10*x[1]/pow(x[0],3));
  g[1]=2*x[0]+(5/pow(x[0],2));
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

  //Problem 1a
  d={1.0,0.0};
  x0={-1.0,0.5};

  double f1=func1a(x0);
  printf("f1=%2.8f\n",f1);

  valarray<double> g1a=grad1a(x0);
  for (int i=0;i<g1a.size();i++) {
      printf("gradient1=%2.8f\n",g1a[i]);
  }

  std::cout<<"Minimize f(x,y)=-12y+4x^2+4y^2+4xy"<<std::endl;
  t=clock();
  std::cout<<"Running Armijio's Rule Inexact Line Search..."<<std::endl;
  ar1=armijo_rule(func1a,grad1a,x0,d,stepsize,eta,theta);
  printf("Final output = [%2.8f %2.8f]\n",ar1[0],ar1[1]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
 
  //Problem 1b
  d={1.0,0.0};
  x0={1.0,1.0};
  std::cout<<"Minimize objective function in problem 1b"<<std::endl;
  t=clock();
  std::cout<<"Running Armijio's Rule Inexact Line Search..."<<std::endl;
  ar2=armijo_rule(func1b,grad1b,x0,d,stepsize,eta,theta);
  printf("Final output = [%2.8f %2.8f]\n",ar2[0],ar2[1]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  
  d={0.0,-1.0};
  x0={3.0,3.0};
  std::cout<<"Minimize objective function in problem 1b"<<std::endl;
  t=clock();
  std::cout<<"Running Armijio's Rule Inexact Line Search..."<<std::endl;
  ar3=armijo_rule(func1b,grad1b,x0,d,stepsize,eta,theta);
  printf("Final output = [%2.8f %2.8f]\n",ar3[0],ar3[1]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

  return 0;
}
