#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void)
{

//  srand48(time(NULL));
  srand48(1.);
  printf("%.12f is a random number in [0.0, 1.0).\n", drand48());

  exit(0);
}


