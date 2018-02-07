/* isfinite example */
#include <stdio.h>      /* printf */
#include <math.h>       /* isfinite, sqrt */
#include <assert.h>     /* assert */
#include <algorithm>    // std::max
#include <iostream>
using namespace std;
int main()
{
 int x = 1.;
 int y = 2.;
 int z = 3.;
 
 int out = std::max(x, std::max(y,z));
 printf ("max(x,max(y,z))      : %d\n",out);

 float tmp1 = 1./0.;
 cout<<"tmp1="<<tmp1<<" Inf*0="<< tmp1 * 0.; 
 printf ("tmp1=%d, tmp1*0 = %d  \n",tmp1, tmp1*0.);
  return 0;
}
