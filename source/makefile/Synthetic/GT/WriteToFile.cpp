#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>     /* assert */
#include <stdlib.h>     /* srand, rand */
#include <math.h>
#include <time.h>       /* time */

int main ()
{
    FILE * pFile;
    float buffer[3] = { 0.34,  0.58,  0.48 }; // center of mass of tumour
    
    std::cout << buffer[0] << " " << buffer[1] << " " << buffer[2] << std::endl;
    
    pFile = fopen ("HGG_TumorIC.bin", "wb");
    fwrite (buffer , sizeof(float), sizeof(buffer), pFile);

    // terminate
    fclose (pFile);
   // free (buffer);


  return 0;
}


