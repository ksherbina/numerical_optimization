//Implementation of Golden Section Search to find a local minimum of a user-defined function
//for both 1-dimensional and multidimensional inputs.
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>
#include <typeinfo>

using std::valarray;

valarray<double> golden_section_search(int n, double (*f)(valarray<double>), valarray<double> x, double a, double b, valarray<double> d, double e)
{
  // User must supply the function f: R^n -> R.
  valarray<double> xs (n), xu (n), xr (n), xl (n);
  double diff, diff0, fxl, fxr;
  // Since you want to stay within the initial left and right endpoints, use the solution to
  //t^2+t-1=0 such that t>0 as the golden ratio r.
  double r=(0.5)*(sqrt(5.0)-1);
  //printf("golden ratio = %2.8f \n",r);
  xs=x+a*d; //initialize left endpoint
  xu=x+b*d; //initialize right endpoint
  xr=xs+r*(b-a)*d; //initialize the right-hand search point
  xl=xs+(1.0-r)*(b-a)*d; //initialize the left-hand search point
  fxr=f(xr);
  fxl=f(xl);
  diff=sqrt(pow(xu-xs, 2).sum());
  
  /*
  int s1=18,s2=16;
  printf("***NOTE: xs=left endpoint; xu=right endpoint; xl=left-hand search point; xr=right-hand search point***\n");
  //printf("%s & %s & %*s & %*s & %*s & %*s & %*s \\\\ \n","Iteration","xs",s2,"xu",s2,"xl",s2,"xr",s1,"f(xl)",s2,"f(xr)");
  //printf("%d & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f \\\\ \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
  printf("%s %s %*s %*s %*s %*s %*s\n","Iteration","xs",s2,"xu",s2,"xl",s2,"xr",s1,"f(xl)",s2,"f(xr)");
  printf("%d %*.8f %*.8f %*.8f %*.8f %*.8f %*.8f \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
  */
  int i=0;
  // invariant: The euclidean norm of the difference between the left and right
  // endpoints is greater than the stopping tolerance
  while (diff>e) {
    diff0=diff;
    i++;
    if (fxr>fxl) {
      //printf("fxr greate than fxl\n");
      xu=xr;
      xr=xl;
      fxr=fxl;
      xl=xs+(xu-xr);
      fxl=f(xl);
    } else {
      //printf("fxr less than or equal to fxl\n");
      xs=xl;
      xl=xr;
      fxl=fxr;
      xr=xu-(xl-xs);
      fxr=f(xr);
    }
    /*
    printf("%d %*.8f %*.8f %*.8f %*.8f %*.8f %*.8f \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
    */
    diff=sqrt(pow(xu-xs, 2).sum());
    if (diff0==diff) {
      printf("Cannot find a minimum given the interval of uncertainty of %2.8f.\n",e);
      printf("||xu-xs|| = %2.8f\n",diff);
      break;
    }
  }

  return (xu+xs)/2;
}
