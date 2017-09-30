//Implementation of Dichotomous Search to find a local minimum of a user-defined function
//for both 1-dimensional and multidimensional inputs.
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>
#include <typeinfo>

using std::valarray;

valarray<double> dichotomous_search(int n, double (*f)(valarray<double>), valarray<double> x, double a, double b, valarray<double> d, double epsilon, double eta)
{
  //User must supply the function f: R^n -> R.
  valarray<double> xs (n), xu (n), xr (n), xl (n);
  double diff, diff0, fxl, fxr;
  //delta is the distinguishability constant such that eta = 2*delta > 0
  printf("Distinguishability constant delta such that %2.8f = 2*delta > 0 :\n",eta);
  valarray<double> delta (eta/2.0, n);
  for (int i=0; i<delta.size(); i++) {
    std::cout<<delta[i]<<' ';
  }
  std::cout<<std::endl;
  
  xs=x+a*d; //initialize left endpoint
  xu=x+b*d; //initialize right endpoint
  xl=0.5*(a+b)*d-delta; //initialize the left-hand search point
  xr=0.5*(a+b)*d+delta; //initialize the right-hand search point
  fxr=f(xr);
  fxl=f(xl);
  diff=sqrt(pow(xu-xs, 2).sum());
  
  int s1=18,s2=16;
  int i=0;
  printf("***NOTE: xs=left endpoint; xu=right endpoint; xl=left-hand search point; xr=right-hand search point***\n");
  //printf("%s & %s & %*s & %*s & %*s & %*s & %*s \\\\ \n","Iteration","xs",s2,"xu",s2,"xl",s2,"xr",s1,"f(xl)",s2,"f(xr)");
  //printf("%d & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f \\\\ \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
  printf("%s %s %*s %*s %*s %*s %*s\n","Iteration","xs",s2,"xu",s2,"xl",s2,"xr",s1,"f(xl)",s2,"f(xr)");
  printf("%d %*.8f %*.8f %*.8f %*.8f %*.8f %*.8f \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
         
  //invariant: The euclidean norm of the difference between the left and right
  //endpoints is greater than the stopping tolerance
  while (diff>epsilon) {
    diff0=diff;
    i++;
    if (fxr>fxl) {
      //printf("fxr greate than fxl\n");
      xu=xr;
    } else {
      //printf("fxr less than or equal to fxl\n");
      xs=xl;
    }
    xl=0.5*(xs+xu)*d-delta;
    xr=0.5*(xs+xu)*d+delta;
    fxl=f(xl);
    fxr=f(xr);
    //printf("%d & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f & %*.8f \\\\ \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
    printf("%d %*.8f %*.8f %*.8f %*.8f %*.8f %*.8f \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
    diff=sqrt(pow(xu-xs, 2).sum());
    if (diff0==diff) {
      printf("Cannot find a minimum given the interval of uncertainty of %2.8f.\n",epsilon);
      printf("||xu-xs|| = %2.8f\n",diff);
      break;
    }
  }
  return (xu+xs)/2;
}
