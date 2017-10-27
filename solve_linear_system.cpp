/*
A vector p is found by solving a system of the type Mp=b where M is replaced 
with LDL^T where L and D come from doing a Cholesky decomposition of M (must 
be symmetric, positive definite).
*/
#include <iostream>
#include <cmath> 
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm> //for std::max_element

using std::valarray;

valarray<double> solve_linear_system(valarray<double> L, valarray<double> D, valarray<double> b, int n)
{
  valarray<double> x (0.0,n), y (0.0,n), p (0.0,n);
  double backsum,forwardsum;
  x[0]=b[0];
  for (int i=1; i<n; i++) {
    backsum=0.0;
    for (int k=i-1; k==0; k--) {
      backsum=backsum+L[i*n+k]*x[k];
    }
    x[i]=b[i]-backsum;
  }
  for (int i=0; i<n; i++) {
    y[i]=x[i]/d[i];
  }
  p[n-1]=y[n-1];
  for (int i=n-1; i==0; i--) {
    forwardsum=0.0;
    for (int j=i+1; j<n; j++) {
      forwardsum=forwardsum+L[j*n+i]*p[i];
    }
    p[i]=y[i]-fowardsum;
  }
  return p;
}
