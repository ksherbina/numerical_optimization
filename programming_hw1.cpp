//Function, inputs, and parameters for problem 2a of homework 1 in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "golden_section_search.h"
#include "dichotomous_search.h"

using std::valarray;

double func1(valarray<double> x) {
    //user-defined function f such that f:R^n->R.
    return (pow(x,3.0)-3.0*(pow(x,2.0)))[0];
}

double func2(valarray<double> x) {
  //user-defined function f such that f:R^n->R.
  valarray<double> z = x;
  for (int i=0;i<x.size();i++) {
    z[i] = 1/x[i];
  }
  return (20*z+pow(x,2.0))[0];
}

int main()
{
  clock_t t;
  int n=1; //User must specify number of dimensions
  valarray<double> d (n), x (n), xs (n), gss (n), ds (n);
  double a,b,error,delta;
  
  //Problem 2a
  d[0]=1.0;
  x[0]=0.0;
  a=0;
  b=3;
  error=pow(10.0,-4); //Interval of Uncertainty
  
  std::cout<<"Minimize f(x)=x^3-3*x^2"<<std::endl;
  t=clock();
  std::cout<<"Running Golden Section Search..."<<std::endl;
  gss=golden_section_search(n,func1,x,a,b,d,error);
  printf("Final output = %2.8f\n",gss[0]);
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
  
  return 0;
}
