//Function for problem 2a of homework 1 in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "golden_section_search.h"

using std::valarray;

double func(valarray<double> x) {
    //user-defined function f such that f:R^n->R.
    return (pow(x,3.0)-3.0*(pow(x,2.0)))[0];
}

int main()
{
  clock_t t;
  t=clock();
  int n=1;
  valarray<double> d (n), x (n), xs (n), soln(n);
  double a,b,error,fval1;

  d[0]=1.0;
  x[0]=0.0;
  a=0;
  b=3;
  error=pow(10.0,-4);
  printf("Interval of Uncertainity = %2.8f \n",error);
  //fval1=func(x+b*d);
  //printf("x = %2.8f \n",(x+b*d)[0]);
  //printf("f(x) = %2.8f \n",fval1);
  soln=golden_section_search(n,func,x,a,b,d,error);
  std::cout<<"size of soln: "<<soln.size()<<std::endl;
  printf("soln: %2.8f \n",soln[0]);
  t=clock()-t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  return 0;
}
