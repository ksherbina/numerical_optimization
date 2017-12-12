#ifndef SOLVE_QP_FROTH_H
#define SOLVE_QP_FROTH_H
 
#include <valarray>
using std::valarray;

struct QpSolutions {
    valarray<double> direction;
    valarray<double> multipliers;
};

QpSolutions solve_qp_froth(double fxk, valarray<double> gxk, valarray<double> hxk, valarray<double> B, valarray<double> c, 
         valarray<double> cg0, valarray<double> cg1, valarray<double> ch0, valarray<double> ch1) {
 
#endif
