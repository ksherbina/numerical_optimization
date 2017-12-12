/*
// min_x max_i \sum_i |f_i(x)| where 
f_i(x) come from the Freudenstein and Roth function.

To use SQP, this problem was converted to 
min z where z = max_i \sum_i |f_i(x)|
subject to f_0 - z_0 <= 0
           f_1 - z_0 <= 0
           -f_0 - z_1 <= 0
           -f_1 - z_1 <= 0
*/


#include <iostream>
#include <valarray>
#include <ilcplex/ilocplex.h>

using std::valarray;

struct QpSolutions {
    valarray<double> direction;
    valarray<double> multipliers;
};


QpSolutions solve_qp_froth(double fxk, valarray<double> gxk, valarray<double> hxk, valarray<double> B, valarray<double> c, 
         valarray<double> cg0, valarray<double> cg1, valarray<double> ch0, valarray<double> ch1) {
    /*
    c is an array of the constraints evaluated at the current point 
    cg0 is the gradient of the first constraint evaluated at the current point
    cg1 is the gradient of the second constraint evaluated at the current point 
    */
    QpSolutions qp;
    IloEnv      env;
    IloModel    mod0 (env);
    
    IloNumVarArray d (env, ILOFLOAT);
    
    p.setNames("d");

    IloExpr     expr (env);
    expr += fxk + gxk[0] * d[0] + gxk[1] * d[1]
    expr += 0.5 * (d[0] * (d[0] * hxk[0] + d[1] * hxk[2]) + d[1] * (d[0] * hxk[1] + d[1] * hxk[3]));
    
    // Constraint 1 & 2
    // Option 1
    mod.add( c[0] + cg0[0] * d[0] + cg0 * d[0] - expr <= 0 );
    mod.add( c[1] + cg1[0] * d[0] + cg1 * d[0] - expr <= 0 );
    mod.add( -c[0] - cg0[0] * d[0] - cg0 * d[0] - expr <= 0 );
    mod.add( c[1] - cg1[0] * d[0] - cg1 * d[0] - expr <= 0 );
    
    
    // Objective Function
    // x1 + 2 x2 + 3 x3
    //            - 0.5 ( 33*x1*x1 + 22*x2*x2 + 11*x3*x3
    //                             - 12*x1*x2 - 23*x2*x3 )
    IloExpr obj (env);
    obj == expr;
    
    mod.add( IloMinimize(env, obj) );
    
    IloCplex cpx (env);
    cpx.extract( mod );
    
    //cpx.setOut( env.getNullStream() );    // no output
    
    cpx.exportModel("QP.lp");
    
    cpx.solve();
    
    cout << cpx.getObjValue() << endl;
    // Get the direction
    for (int i=0; i<gxk.size(); i++) {
      qp.direction[i] = cpx.getValue( d[i] );
    }
    // Get the multipliers
    status = CPXXgetpi (env, cpx, pi, 0, CPXXgetnumrows(env,cpx)-1);
    for (int i=0; i<gxk.size(); i++) {
      qp.multipliers[i] = pi[i];
    }
    
    cout << cpx.getValue( expr ) << endl;
    
    env.end();
    
    return qp;
}
