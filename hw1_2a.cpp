//Function for problem 2a of homework 1 in ISE 520 Fall 2017
#include <iostream>
#include <stdio.h>
#include <time.h>  //for clock_t, clock, CLOCKS_PER_SEC

int main()
{
  clock_t t;
  t = clock();
  int n = 1;
  std::cout<<n<<std::endl;
  float d[n] = {1};

  float a = 0.5;
  float b = 1.5;
  printf ("%2.4f %2.4f \n",a,b);
  float x0[n] = {2};
  t = clock() - t;
  printf ("Runtime of algorithm: %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
  return 0;
}
