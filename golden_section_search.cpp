//Implementation of Golden Section Search in order to find minimum of a 
//quasi-convex function.
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>

using std::valarray;

valarray<double> golden_section_search(int n, valarray<double> x, double a, valarray<double> d)
{
  valarray<double> xs (n);
  xs = x + a*d;
  printf("left endpoint: %2.8f \n",xs[0]);
  return xs;
}
