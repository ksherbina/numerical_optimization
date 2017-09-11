//Implementation of Golden Section Search in order to find minimum of a 
//quasi-convex function.
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>

using std::valarray;

valarray<double> golden_section_search(int n, valarray<double> (*f)(valarray<double>), valarray<double> x, double a, double b, valarray<double> d, double e)
{
  valarray<double> xs (n), xu (n), xr (n), xl (n), fxr (n), fxl (n);
  double diff;
  double r=(0.5)*(1.0+sqrt (5.0));
  printf("golden ratio = %2.8f \n",r); 
  xs=x+a*d; //initialize left endpoint
  xu=x+b*d; //initialize right endpoint
  std::cout<<"size of xs: "<<xs.size()<<" size of xu: "<<xu.size()<<std::endl;
  printf("xs: %2.8f, xu: %2.8f \n",xs[0],xu[0]);
  xr=xs+r*(b-a)*d;
  xl=xs+(1.0-r)*(b-a)*d;
  fxr=f(xr);
  fxl=f(xl);
  std::cout<<"size of f(xr): "<<fxr.size()<<" size of f(xl): "<<fxl.size()<<std::endl;
  printf("xr = %2.8f, f(xr) = %2.8f \n",xr[0],fxr[0]);
  printf("xl = %2.8f, f(xl) = %2.8f \n",xl[0],fxl[0]);
  diff=sqrt(pow (xu-xs, 2).sum());  
  printf("invariant: %2.8f \n",diff);

  return xs;
}
