//Function, inputs, and parameters for problem 1 of homework 2a in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
//#include "armijo_rule.h"

using std::valarray;

double func2a(double x) {
    //user-defined function f such that f:R->R.
    return pow(x,3)-3*pow(x,2);
}

/*
valarray<double> grad1a(valarray<double> x) {
    //gradient of func1a; in R^n
    valarray<double> g (2);
    g[0]=8.0*x[0]+4.0*x[1];
    g[1]=4.0*x[0]+8.0*x[1]-12.0;
    return g;
}
*/

int main()
{
  clock_t t;
  int n=2; //User must specify number of dimensions
  
  double x0,a,b;
  double epsilon=pow(10.0,-8);
  double theta=0.5;

  //Problem 1a
  x0=1.2;
  a=1.0;
  b=4.0;

  double f1=func2a(x0);
  printf("f1=%2.8f\n",f1);
  /*
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
