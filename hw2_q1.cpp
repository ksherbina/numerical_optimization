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
/*
double func2(valarray<double> x) {
  //user-defined function f such that f:R^n->R.
  valarray<double> z = x;
  for (int i=0;i<x.size();i++) {
    z[i] = 1/x[i];
  }
  return (20*z+pow(x,2.0))[0];
}
*/
int main()
{
  clock_t t;
  int n=2; //User must specify number of dimensions

  valarray<double> d(n),x0(n),ar1(n),ar2(n);
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
  ar1=armijo_rule(n,func1a,g1a,x0,d,stepsize,eta,theta);
  printf("Final output = [%2.8f %2.8f]\n",ar1[0],ar1[1]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
 
  /*
  //User specified parameters to run dichotomous search different from golden section search
  error=2*pow(10.0,-2); //Interval of Uncertainty
  delta=pow(10.0,-2); //Distinguishability constant for dichotomous search
  
  t=clock();
  std::cout<<std::endl<<"Running Dichotomous Search..."<<std::endl;
  ds=dichotomous_search(n,func1,x,a,b,d,error,delta);
  printf("Final output = %2.8f\n",ds[0]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

  //Problem 2b
  x[0]=0.01;
  a=0.04;
  b=6;
  error=pow(10.0,-4); //Interval of Uncertainty
  std::cout<<std::endl<<"Minimize f(x)=(20/x)+x^2"<<std::endl;
  t=clock();
  std::cout<<"Running Golden Section Search..."<<std::endl;
  gss=golden_section_search(n,func2,x,a,b,d,error);
  printf("Final output = %2.8f\n",gss[0]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  
  error=2*pow(10.0,-2); //Interval of Uncertainty
  delta=pow(10.0,-2); //Distinguishability constant for dichotomous search
  
  t=clock();
  std::cout<<std::endl<<"Running Dichotomous Search..."<<std::endl;
  ds=dichotomous_search(n,func2,x,a,b,d,error,delta);
  printf("Final output = %2.8f\n",ds[0]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
*/
  return 0;
}
