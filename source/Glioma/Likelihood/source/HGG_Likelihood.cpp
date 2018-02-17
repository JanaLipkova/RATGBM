//
//  HGG_Likelihood.cpp
//  GliomaXcode
//
//  Created by Lipkova on 15/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#include "HGG_Likelihood.h"

HGG_Likelihood::HGG_Likelihood(const int argc, const char ** argv):parser(argc, argv)
{
    ifstream mydata("LikelihoodInput.txt");
    
    if (mydata.is_open())
    {
        mydata >> T1uc;
        mydata >> T2uc;
        mydata >> slope;
        mydata.close();
    }
    else
    {
        printf("Aborting: missing input file LikelihoodInput.txt \n");
        abort();
    }
    
    printf("MRI: T1uc=%f, T2uc =%f, slope=%f \n", T1uc, T2uc, slope);
}


long double HGG_Likelihood::_computeTiLogLikelihood(MatrixD3D model, int day, int Ti)
{
    /* 1) Read in Ti data */
    char filename[256];
    sprintf(filename,"T%dw_D%02d.dat",Ti,day);

    
    MatrixD3D data(filename);
    int dataX = data.getSizeX();
    int dataY = data.getSizeY();
    int dataZ = data.getSizeZ();
    
    MatrixD3D mask("Mask.dat");
    int Npoints = 0;
    
    long int N = dataX*dataY*dataZ;
    assert(N == model.getSizeX() * model.getSizeY() * model.getSizeZ() );
    
    long double sum = 0.;
    int cor_leng = 6;
    
    for (int iz = 0; iz < dataZ; iz=iz+cor_leng )
        for (int iy = 0; iy < dataY; iy++)
            for (int ix = 0; ix < dataX; ix++)
            {
                if(mask(ix,iy,iz) > 0.01)
                {
                    sum += _computeLogBernoulli(model(ix,iy,iz), data(ix,iy,iz), Ti);
                    Npoints++;
                }
            }
    
    printf("LogLike of day %d of T%d = %Lf (using Npoints=%d) \n", day, Ti, sum, Npoints);
    return sum;
}


/* Likelihood based on Bernoulli distr. of double logsitic sigmoid
 1) compute alpha
 - should be in (0,1)
 - rounding errros can make alpha = 0 or =1
 - if that is the case, correct it since it will be used in log() */
long double HGG_Likelihood::_computeLogBernoulli(double u, double y, int Ti)
{
    double uc, is2;
    
    if(Ti == 1){
        uc = T1uc;
        is2 = 1./slope;
    }
    else{
        uc = T2uc;
        is2 = 1./slope;
    }
    
    double diff = u - uc;
    
    long double alpha = 0.5 + 0.5 * sgn(diff) * (1. - exp( -diff * diff * is2));
    return  (y == 1 ) ? log(alpha) : log(1.-alpha);
}

int HGG_Likelihood::sgn(double d)
{
    double eps = 0.0;
    if (d < -eps) { return -1; }
    else { return d > eps; }
}

void HGG_Likelihood::_writeToFile(long double output)
{
    long double MinusLogLikelihood = - output;
    
    FILE* myfile = fopen("Likelihood.txt", "w");
    if (myfile!=NULL)
        fprintf(myfile, "%Lf \n", MinusLogLikelihood);
    
    fclose(myfile);
}

void HGG_Likelihood::run()
{
    char filename[256];
    sprintf(filename,"M_UQ_J09.dat");
    MatrixD3D model_J9(filename);
    
    sprintf(filename,"M_UQ_J1.dat");
    MatrixD3D model_J11(filename);

    
    int day = 9;
    long double LT1_J9  = _computeTiLogLikelihood(model_J9, day, 1);
    long double LT2_J9  = _computeTiLogLikelihood(model_J9, day, 2);
    
    day = 11;
    long double LT1_J11  = _computeTiLogLikelihood(model_J11, day, 1);
    long double LT2_J11  = _computeTiLogLikelihood(model_J11, day, 2);
    
    long double costFunction = LT1_J9 + LT2_J9 + LT1_J11 + LT2_J11 ;
   
    printf("LT1_J9=%Lf, LT2_J9=%Lf, LT1_J11=%Lf, LT2_J11=%Lf, \n", LT1_J9, LT2_J9, LT1_J11, LT2_J11);
    printf("LogLike = %Lf \n", costFunction);
    _writeToFile(costFunction);

}
