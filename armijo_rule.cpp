/*
Implementation of Inexact Search using Armijo's rule (refer to Chapter 8 in
Bazaraa, M.S., Sherali, H.D., & Shetty, C.M. "Nonlinear Programming (3rd ed)" 
to minimize a function.
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

int comparison (int iter,double (*f)(valarray<double>),const double& grad,const valarray<double>& xc,const valarray<double>& p, double multiplier, double th)
{
  double fxc,lc,fxk;
  valarray<double> xk;
  fxc=f(xc);
  lc=fxc+multiplier*th*grad;
  xk=xc+(multiplier*p);
  fxk=f(xk);

  int s1=16;
  //printf("%d & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f \\\\ \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
  //printf("%d & [%.8f %.8f] & [%.8f %.8f] & %.8f & %.8f \\\\ \n",iter,xc[0],xc[1],xk[0],xk[1],fxc,lc);
  printf("%d %.8f %*.8f %*.8f\n",iter,multiplier,s1,fxk,s1,lc);
  
  if (fxk<=lc) {
    return 1;
  } else {
    return 0;
  }
}

double armijo_rule(double (*f)(valarray<double>),valarray<double> (*g)(valarray<double>),valarray<double> x0,valarray<double> d, double eta, double theta)
{
  //User must supply the function f: R^n -> R and the gradient of the function g:R^n.
  double alpha,step,descent;
  int t=1;
  valarray<double> gx0;
  
  alpha=1.0;
  gx0=g(x0);
  descent=std::inner_product(std::begin(gx0),std::end(gx0),std::begin(d),0.0);
  //printf("Initial descent direction is %.8f",descent);
  if (descent<0) {
    printf("NOTE: LH is the left hand side of Armijo's rule and RH is the right hand side\n");
    printf("%s %s %s %s\n","Iteration","Steplength","LH","RH");
    if (comparison(t,f,descent,x0,d,alpha,theta)==1) {
      printf("Double stepsize");
      step=pow(eta,t)*alpha;
      while (comparison(t,f,descent,x0,d,step,theta)==0) {
        alpha=step;
        step=pow(eta,t)*alpha;
        t+=1;
      }
    } else {
      //std::cout<<"f(x_i) > r(multiplier)"<<std::endl;
      printf("Halve stepsize");
      step=pow(1/eta,t)*alpha;
      printf("Current step size=%.8f",step);
      while (comparison(t,f,descent,x0,d,step,theta)==0) {
        alpha=step;
        step=pow(1/eta,t)*alpha;
        t+=1;
        std::cout<<step<<std::endl;
      }
    }
    printf("Acceptable steplength: %.8f",step);
    return step;
  } else {
    printf("Current direction is not a descent direction");
    return 0.0;
  }
}
