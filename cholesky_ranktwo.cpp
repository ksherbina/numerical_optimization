/*
Find Cholesky factors of a square matrix by first adding 
and subtracting a rank one matrix to current matrix.
*/
#include <limits>
#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm> //for std::max_element
#include "cholesky.h" //for std::max_element
#include "solve_linear_system.h"

using std::valarray;

CholeskyFactors cholesky_ranktwo(valarray<double> v, valarray<double>u, CholeskyFactors chol0) {
  CholeskyFactors result;
  int n = v.size();
  valarray<double> d(0.0, n), l(0.0, n*n), vmat(0.0, n*n), pvec(0.0, n);
  double p, b, backsum;
  for (int i = 0; i < n; i++) {
    vmat[i] = v[i];
  }
  double t0 = 1.0;
  double t;
  for (int j = 0; j < n; j++) {
    p = vmat[j * n + j];
    t = t0 + (p / chol0.diagonal[j]);
    d[j] = chol0.diagonal[j] * t / t0;
    b = p / (chol0.diagonal[j] * t);    
    t0 = t;
    for (int r = j + 1; r < n; r++) {
      vmat[r * n + (j + 1)] = vmat[r * n + j] - (p * chol0.lower_triangular[r * n + j]);
      l[r * n + j] = chol0.lower_triangular[r * n + j] + (b * vmat[r * n + (j + 1)]);
    }
  }

  valarray<double> dn(0.0, n), ln(0.0, n*n), umat(0.0, n*n);
  for (int i = 0; i < n; i++) {
    umat[i] = u[i];
  }

  pvec[0] = u[0];
  for (int i = 1; i < n; i++) {
    backsum = 0.0;
    for (int k = i - 1; k >= 0; k--) {
      backsum = backsum + l[i * n + k] * pvec[k];
    }
    pvec[i] = u[i] - backsum;
  }

  double pdp = 0.0;
  for (int i = 0; i < n; i++) {
    pdp += pow(pvec[i], 2)*d[i];
  }

  double tn = 1 - pdp;

  if (tn < std::numeric_limits<double>::epsilon()) {
    tn = std::numeric_limits<double>::epsilon();
  }
  
  for (int j = n - 1; j < 0; j--) {
    t = tn + (pow(pvec[j], 2) / d[j]);
    dn[j] = d[j] * tn / t;
    b = -pvec[j] / (d[j] * tn);
    umat[j * n + j] = pvec[j];
    tn = t;
    for (int r = j + 1; r < n; r++) {
      ln[r * n + j] = l[r * n + j] + (b * umat[r * n + (j + 1)]);
      umat[r * n + j] = umat[r * n + (j + 1)] + (pvec[j] * l[r * n + j]);
    }
  }
  result.lower_triangular = ln;
  result.diagonal = dn;

  return result;
}
