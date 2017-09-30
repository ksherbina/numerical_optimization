//Implementation of Dichotomous Search to find a local minimum of a user-defined function
//for both 1-dimensional and multidimensional inputs.
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>
#include <tuple>

using std::valarray;

std::tuple<double, double>  newton_1d(double (*f)(double),double (*g)(double),double (*h)(double),double x0, double a, double b, double epsilon, double theta)
{
  //User must supply the functions f: R^n -> R, the gradient of f g:R^n, and hessian of f h:R^n.
  double d,alpha,xd,fxd;
  int iter=1;
  double xc=x0;
  double fxc=f(xc);
  double gxc=g(xc);
  double hxc=h(xc);
  

  int s1=22,s2=16,s3=14;
  //printf("%s & %s & %*s & %*s & %*s & %*s & %*s \\\\ \n","Iteration","xs",s2,"xu",s2,"xl",s2,"xr",s1,"f(xl)",s2,"f(xr)");
  //printf("%d & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f \\\\ \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
  printf("%s %s %*s %*s %*s\n","Iteration (i)","x_i",s2,"f(x_i)",s2,"g(x_i)",s3,"h(x_i)");
  printf("%d %*.8f %*.8f %*.8f %*.8f\n",iter,s1,xc,s3,fxc,s3,gxc,s3,hxc);
  
  while (std::abs(gxc)>epsilon) {
    if (hxc>0.0) {
      d=-gxc/hxc;
    } else {
      d=-gxc;
    }
    if (d<0.0) {
      alpha=std::min(1.0,(a-xc)/d);
    }
    if (d>0.0) {
      alpha=std::min(1.0,(b-xc)/d);
    }
    xd=xc+alpha*d;
    fxd=f(xd);
    while (fxd>=fxc) {
      alpha=theta*alpha;
      xd=xc+alpha*d;
      fxd=f(xd);
    }
    xc=xd;
    fxc=fxd;
    gxc=g(xc);
    hxc=h(xc);
    iter+=1;
    printf("%d %*.8f %*.8f %*.8f %*.8f\n",iter,s1,xc,s3,fxc,s3,gxc,s3,hxc);
  }
  return std::make_tuple(xc,fxc);
}
