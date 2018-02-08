/* sin example */
#include <stdio.h>      /* printf */
#include <math.h>       /* sin */
#define PI 3.14159265

int main ()
{
  double ICx = 0.315; 
  double ICy = 0.67;
  double ICz = 0.5;

  double CMx = 0.317134;
  double CMy = 0.667676;
  double CMz = 0.494630;

  double x = ICx - CMx;
  double y = ICy - CMy;
  double z = ICz - CMz;

  double r      = sqrt( x*x + y*y + z*z);
  double alpha  = atan (y/x) ;
  double beta   = acos (z/r) ;
  
  printf("In radians: \n r=%f , \n alpha=%f, \n beta=%f, \n", r, alpha, beta );
  
   alpha = alpha * 180. /PI;
   beta = beta * 180. /PI;
 
  printf("In degrees: \n r=%f , \n alpha=%f, \n beta=%f, \n", r, alpha, beta );
 
  printf( "mod(7.1, 2) = %f \n", fmod(7.1, 2.) );
  
  alpha = fmod( (alpha + 360.), 360.);
  beta =  fmod( (beta + 360.), 360.);

   printf("In degrees : \n r=%f , \n alpha=%f, \n beta=%f, \n", r, alpha, beta );
  return 0;
}
