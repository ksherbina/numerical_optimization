#include <iostream>
#include <cmath>
#include <stdio.h> //for printf
#include <valarray>
#include <numeric> //for std::inner_product
#include <string>
#include <iostream>
#include "cholesky.h"
#include "newton.h"
#include "cholesky_ranktwo.h"
#include "solve_linear_system.h"
#include "armijo_rule.h"

using std::valarray;
using std::string;

double scalar_product_of_2_vectors(valarray<double> u, valarray<double> v) {
  /*
  Multiply a 1xn vector u by a nx1 vector v 
  */
  double sum = 0.0;
  for (int j = 0; j < u.size(); j++) {
    sum += (u[j] * v[j]);
  }
  return sum;
}

valarray<double> multiply_matrix_by_vector(valarray<double> M, valarray<double> z) {
  /*
  Multiply a nxn matrix M by a nx1 vector z 
  */
  valarray<double> prod(0.0, z.size());
  double sum;
  for (int j = 0; j < M.size(); j += z.size()) {
    sum = 0.0;
    for (int k = 0; k < z.size(); k++) {
      sum += (M[j + k] * z[k]);
    }
    prod[j/z.size()] = sum;
  }
  return prod;
}

valarray<double> multiply_vector_by_its_transpose(valarray<double> y) {
  valarray<double> mat(0.0, y.size() * y.size());
  for (int j = 0; j < y.size(); j++) {
    for (int k = 0; k < y.size(); k++) {
      mat[j*y.size()+k] = y[j] * y[k];
    }
  }
  return mat;
}

valarray<double> damped_bfgs(valarray<double> B0, valarray<double> delta, valarray<double> gamma) { 

  double theta;
  valarray<double> rk, B;
  
  // compute theta as in equation (18.15) in Nocedal & Wright
  // compute rk as in Procedure 18.2 in Nocedal & Wright 
   
  add_rank1 = multiply_vector_by_its_transpose(r) / scalar_product_of_2_vectors(delta, r);
  hdelta = multiply_matrix_by_vector(B0, delta);
  subtract_rank1 = hdelta / scalar_product_of_2_vectors(delta, hdelta);
  B = B0 - subtract_rank1 + add_rank1

  return B;
}
