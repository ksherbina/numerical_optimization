//Function, inputs, and parameters for problem 2a of homework 1 in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "golden_section_search.h"
#include "dichotomous_search.h"

using std::valarray;

double func(valarray<double> x) {
    //user-defined function f such that f:R^n->R.
    return (pow(x,3.0)-3.0*(pow(x,2.0)))[0];
}

int main()
{
  clock_t t;
  int n=1; //User must specify number of dimensions
  valarray<double> d (n), x (n), xs (n), gss (n), ds (n);
  double a,b,error,delta;
  
  d[0]=1.0;
  x[0]=0.0;
  a=0;
  b=3;
  error=pow(10.0,-4); //Interval of Uncertainty
  delta=pow(10.0,-2); //Distinguishability constant for dichotomous search
  
  t=clock();
  std::cout<<"Running Golden Section Search..."<<std::endl;
  gss=golden_section_search(n,func,x,a,b,d,error);
  printf("Local minimum = %2.8f (given interval of uncertainty of %2.8f)\n",gss[0],error);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  
  //User specified inputs and parameters to run dichotomous search:
  t=clock();
  std::cout<<std::endl<<"Running Dichotomous Search..."<<std::endl;
  ds=dichotomous_search(n,func,x,a,b,d,error,delta);
  printf("Local minimum = %2.8f (given interval of uncertainty of %2.8f)\n",ds[0],error);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  
  return 0;
}
