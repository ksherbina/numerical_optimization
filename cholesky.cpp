/*
Find Cholesky factors of a square matrix and check if the 
matrix is positive definite.
*/
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>

using std::valarray;

double cholesky(valarray<double> A, int n)
{
  valarray<double> C (0.0,n*n), D (0.0,n*n), L (0.0,n*n);
  double sum1,sum2;
  for (int k=0; k<n; k++) {
    L[k*n+k]=1.0;
  }
  D[0]=A[0];
  for (int k=1; k<n; k++) {
    L[k*n]=A[k*n]/D[0];
  }
  for (int j=1; j<n; j++) {
    sum1=0.0;
    for (int s=0; s<j; s++) {
      sum1=sum1+D[s]*(L[j*n+s]*L[j*n+s]);
    } 
    C[j*n+j]=A[j*n+j]-sum1;
    D[j]=C[j*n+j];
    for (int i=j+1; i<n; i++) {
      sum2=0.0;
      for (int t=0; t<j; t++) {
        sum2=sum2+D[t]*L[i*n+t]*L[j*n+t];
      }
      C[i*n+j]=A[i*n+j]-sum2;
      L[i*n+j]=C[i*n+j]/D[j];
    }
  }
  for (int i=0;i<n*n;i++) {
      printf("lower triangular=%2.8f\n",L[i]);
  }
  for (int i=0;i<n;i++) {
      printf("diagonal=%2.8f\n",D[i]);
  }
  return 0.0;
}
