/* isfinite example */
#include <stdio.h>      /* printf */
#include <math.h>       /* isfinite, sqrt */
#include <assert.h>     /* assert */

int main()
{
  printf ("isfinite(0.0)       : %d\n",isfinite(0.0));
  printf ("isfinite(1.0/0.0)   : %d\n",isfinite(1.0/0.0));
  printf ("isfinite(-1.0/0.0)  : %d\n",isfinite(-1.0/0.0));
  printf ("isfinite(sqrt(-1.0)): %d\n",isfinite(sqrt(-1.0)));
  printf ("!isfinite(0)        : %d\n", !isfinite(0.));
  printf ("!isfinite(1./0.)    : %d\n", !isfinite(1./0.));
  assert( isfinite(1.) == 1. ) ;
  return 0;
}
