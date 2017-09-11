//Function for problem 2a of homework 1 in ISE 520 Fall 2017
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC
#include <valarray>
#include "golden_section_search.h"

using std::valarray;

valarray<double> func(valarray<double> x) {
    return pow (x, 3.0)- 3.0*(pow (x, 2.0));
}

int main()
{
  clock_t t;
  t = clock();
  int n = 1;
  valarray<double> d (n);
  valarray<double> x (n);
  valarray<double> xs (n);
  valarray<double> soln (n);

  double a = 0.5;
  double b = 1.5;
  printf ("%2.4f %2.4f \n",a,b);

  d[0] = 1.0;
  x[0] = 2.0;
  soln = golden_section_search(n,x,a,d);
  std::cout<<"size of soln: "<<soln.size()<<std::endl;
  printf("soln: %2.8f \n",soln[0]);
  t = clock() - t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  return 0;
}
