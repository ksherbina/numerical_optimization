//Implementation of Golden Section Search that can be used for multidimensional cases
//to find a local minimum of a user-defined function.
#include <iostream>
#include <cmath>
#include <stdio.h> //include printf
#include <valarray>
#include <typeinfo>

using std::valarray;
using std::greater;

valarray<double> golden_section_search(int n, double (*f)(valarray<double>), valarray<double> x, double a, double b, valarray<double> d, double e)
{
  //user-defined function f: R^n -> R.
  valarray<double> xs (n), xu (n), xr (n), xl (n);
  double diff, fxl, fxr;
  double r=(0.5)*(sqrt(5.0)-1);
  printf("golden ratio = %2.8f \n",r); 
  xs=x+a*d; //initialize left endpoint
  xu=x+b*d; //initialize right endpoint
  std::cout<<"size of xs: "<<xs.size()<<" size of xu: "<<xu.size()<<std::endl;
  printf("xs: %2.8f, xu: %2.8f \n",xs[0],xu[0]);
  xr=xs+r*(b-a)*d;
  xl=xs+(1.0-r)*(b-a)*d;
  fxr=f(xr);
  fxl=f(xl);
  printf("xr = %2.8f, f(xr) = %2.8f \n",xr[0],fxr);
  printf("xl = %2.8f, f(xl) = %2.8f \n",xl[0],fxl);
  diff=sqrt(pow(xu-xs, 2).sum());
  printf("invariant: %2.8f \n",diff);
  
  int s1=18,s2=15;
  int i=0;
  printf("%s %s %*s %*s %*s %*s %*s\n","Iteration","xs",s2,"xu",s2,"xl",s2,"xr",s2,"f(xl)",s2,"f(xr)");
  printf("%d %*.8f %*.8f %*.8f %*.8f %*.8f %*.8f \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
  //invariant: The euclidean norm of the difference between the left and right
  //endpoints is greater than the stopping tolerance
  while (diff>e) {
    if (fxr>fxl) {
      printf("fxr greate than fxl\n");
      xu=xr;
      xr=xl;
      fxr=fxl;
      xl=xs+(xu-xr);
      fxl=f(xl);
    } else {
      printf("fxr less than or equal to fxl\n");
      xs=xl;
      xl=xr;
      fxl=fxr;
      xr=xu-(xl-xs);
      fxr=f(xr);
    }
    printf("%d %*.8f %*.8f %*.8f %*.8f %*.8f %*.8f \n",i,s1,xs[0],s2,xu[0],s2,xl[0],s2,xr[0],s2,fxl,s2,fxr);
    diff=sqrt(pow(xu-xs, 2).sum());
    printf("invariant: %2.8f \n",diff);
    i++;
  }
  return (xu+xs)/2;
}
