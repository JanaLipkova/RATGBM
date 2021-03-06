//
//  Glioma_RAT_UQ.cpp
//  
//  Created by Lipkova on 08/02/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#include "Glioma_RAT_UQ.h"

static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_RAT_UQ::Glioma_RAT_UQ(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////                 RAT GLIOMA UQ                ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    bAdaptivity = parser("-adaptive").asBool();
    pID         = parser("-pID").asInt();
    ICtype      = parser("-ICtype").asInt(0.);
    
    
    ifstream mydata("HGG_InputParameters.txt");
    Dg, Dw, rho, scale;
    
    if (mydata.is_open())
    {
        mydata >> Dw;
        mydata >> rho;
        mydata >> scale;
        mydata.close();
    }
    
    
    switch (ICtype) {
        case 0:
        {
           if(bVerbose) printf("Initialising point tumour \n");
            _ic_rat_point_tumor(*grid, pID);
        }
            break;
            
        case 1:
        {
            if(bVerbose) printf("Reading tumour from file tumour \n");
	    _ic_anatomy(*grid, pID);
            _readInTumourFromFile(*grid, pID, scale);
        }
            break;
            
        case 2:
        {
            if(bVerbose) printf("Initialising two point tumours \n");
            _ic_rat_two_foci_tumor(*grid, pID);
        }
            break;
            
        case 3:
        {
            _ic_rat_elongated_tumor(*grid, pID);
        }
            break;
            
            
        default:
            break;
    }

    
    if(parser("-bDumpIC").asBool(1))
    _dump(0);
    
    isDone              = false;
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= parser("-dumpstart").asDouble();
    whenToRefine        = parser("refinefreq").asDouble();
    whenToRefineOffset  = whenToRefine;

//    whenToWrite			= whenToWriteOffset;
    numberOfIterations	= 0;
    
}

Glioma_RAT_UQ::~Glioma_RAT_UQ()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
/* Initisalise tumour as a point sources, i.e. small smooth sphere
   1) read in anatomies - rescaled to [0,1]^3
   2) read in tumor center of mass + initialize tumor around
   3) set length of brain */
void Glioma_RAT_UQ::_ic_rat_point_tumor(Grid<W,B>& grid, int pID)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
    sprintf(dataFolder,"/home/baldesi/Glioma/RATGBM/source/Anatomy/F98/");
#endif
    
    sprintf(patientFolder, "%sM%02d/M%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s \n", patientFolder);
    
    sprintf(anatomy, "%s_gm.dat", patientFolder);
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%s_wm.dat", patientFolder);
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%s_csf.dat", patientFolder);
    MatrixD3D CSF(anatomy);
    sprintf(anatomy, "%s_mask.dat", patientFolder);
    MatrixD3D MASK(anatomy);

    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.117;   // voxel spacing 117 µm, convert to mm -> L ~ 14 mm
    
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    std::cout<<"L="<<L<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(_DIM);
    _readInTumorPosition(tumor_ic);
    
    // Tumour initial injection for F98: 3µL = 3 mm^3 --> volume corresponding to sphere with radius r=0.8947 mm
    
    const Real tumorRadius = 0.01;//0.8947/L ; // map to [0,1]^3 space
    const Real smooth_sup  = 2;		// 2.suppor of smoothening, over how many gp to smooth
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const Real h = vInfo[0].h[0];  // make sure inference and propagation is run with same h for IC
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)round( x[0] / brainHx  );
                    int mappedBrainY = (int)round( x[1] / brainHy  );
                    int mappedBrainZ = (int)round( x[2] / brainHz  );
                    
                    
                    Real PGt, PWt, Pcsf, Pmask;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PGt   = 0.;
                        PWt   = 0.;
                        Pcsf  = 0.;
                        Pmask = 0.;
                    }
                    else
                    {
                        PGt  =  GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PWt  =  WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pcsf =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pmask = MASK(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    double all = PGt + PWt + Pcsf;
                    
                    if( all >  0.1 )
                    {
                        
                        if(Pcsf > 0.1)
                        {
                            block(ix,iy,iz).p_csf = 1.0;
                            block(ix,iy,iz).p_g   = 0.0;
                            block(ix,iy,iz).p_w   = 0.0;
                        }
                        else
                        {
                            if(PWt > 0.5)
                            {
                                block(ix,iy,iz).p_w = 1.0;
                                block(ix,iy,iz).p_g = 0.0;
                            }
                            else
                            {
                                block(ix,iy,iz).p_w = 0.0;
                                block(ix,iy,iz).p_g = 1.0;
                            }
                            
                        }
                        
                    }
                    
                    // fill the holes in the anatomy segmentations for the rats
                    all = block(ix,iy,iz).p_csf + block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    
                    if( (Pmask > 0.1 )&&( all< 0.1 )  )
                        block(ix,iy,iz).p_g = 1.;
                    
                    /* tumor */
                    const Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                    const Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    const Real psi = (dist - tumorRadius)*iw;
                    
                    if ((psi < -1)&& ((PGt>0.001) || (PWt >0.001)) )		// we are in tumor
                        block(ix,iy,iz).phi = 1.0;
                    else if(( (-1 <= psi) && (psi <= 1) )&& ((PGt>0) || (PWt >0)) )
                        block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                    else
                        block(ix,iy,iz).phi = 0.0;
                    
                    block(ix,iy,iz).phi = min(1.0, 4.*block(ix,iy,iz).phi ); //to ensure that max concentraiton is 1, for correct front propagation speed
                    block(ix,iy,iz).chi = Pmask;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}

void Glioma_RAT_UQ::_readInTumorPosition(vector<Real>& tumorIC )
{
    typedef float dataType;
    
    FILE* fp;
    fp = fopen("HGG_TumorIC.bin", "rb");
    if (fp == NULL) {fputs ("File error", stderr); exit (1);}
    
    // obtain file size
    fseek (fp , 0 , SEEK_END);
    long int size = ftell (fp);
    rewind (fp);
    
    // allocate memory to contain the whole file:
    dataType * buffer;
    buffer = (dataType*) malloc (sizeof(dataType)*size);
    if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
    
    // copy the file into the buffer:
    size_t result;
    result = fread (buffer, 1, size, fp);
    if (result != size) {fputs ("Reading error",stderr); exit (3);}
    
    for (int i = 0; i < _DIM; ++i)
        tumorIC[i] = (Real)buffer[i];
    
    
    free(buffer);
    fclose (fp);
}



void Glioma_RAT_UQ::_ic_anatomy(Grid<W,B>& grid, int pID)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
    sprintf(dataFolder,"/home/baldesi/Glioma/RATGBM/source/Anatomy/F98/");
#endif
    
    sprintf(patientFolder, "%sM%02d/M%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s \n", patientFolder);
    
    sprintf(anatomy, "%s_gm.dat", patientFolder);
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%s_wm.dat", patientFolder);
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%s_csf.dat", patientFolder);
    MatrixD3D CSF(anatomy);
    sprintf(anatomy, "%s_mask.dat", patientFolder);
    MatrixD3D MASK(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.117;   // voxel spacing 117 µm, convert to mm -> L ~ 14 mm
    
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    std::cout<<"L="<<L<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    

    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)round( x[0] / brainHx  );
                    int mappedBrainY = (int)round( x[1] / brainHy  );
                    int mappedBrainZ = (int)round( x[2] / brainHz  );
                    
                    
                    Real PGt, PWt, Pcsf, Pmask;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PGt   = 0.;
                        PWt   = 0.;
                        Pcsf  = 0.;
                        Pmask = 0.;
                    }
                    else
                    {
                        PGt  =  GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PWt  =  WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pcsf =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pmask = MASK(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    double all = PGt + PWt + Pcsf;
                    
                    if( all >  0.1 )
                    {
                        if(Pcsf > 0.1)
                        {
                            block(ix,iy,iz).p_csf = 1.0;
                            block(ix,iy,iz).p_g   = 0.0;
                            block(ix,iy,iz).p_w   = 0.0;
                        }
                        else
                        {
                            if(PWt > 0.5)
                            {
                                block(ix,iy,iz).p_w = 1.0;
                                block(ix,iy,iz).p_g = 0.0;
                            }
                            else
                            {
                                block(ix,iy,iz).p_w = 0.0;
                                block(ix,iy,iz).p_g = 1.0;
                            }
                            
                        }
                        
                    }
                    
                    // fill the holes in the anatomy segmentations for the rats
                    all = block(ix,iy,iz).p_csf + block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    
                    if( (Pmask > 0.1 )&&( all< 0.1 )  )
                        block(ix,iy,iz).p_g = 1.;
                    
                    block(ix,iy,iz).chi = Pmask;
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}




void Glioma_RAT_UQ:: _readInTumourFromFile(Grid<W,B>& grid, int pID, Real scale)
{
    char dataFolder   [200];
    char patientFolder[200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
   sprintf(dataFolder,"/home/baldesi/Glioma/RATGBM/source/Anatomy/F98/");
#endif
    
    sprintf(patientFolder, "%sM%02d/M%02d_TumourIC.dat",dataFolder,pID,pID);
    printf("Reading tumour from %s \n", patientFolder);
    MatrixD3D Tumor(patientFolder);
    
    
    int brainSizeX = (int) Tumor.getSizeX();
    int brainSizeY = (int) Tumor.getSizeY();
    int brainSizeZ = (int) Tumor.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)round( x[0] / brainHx  );
                    int mappedBrainY = (int)round( x[1] / brainHy  );
                    int mappedBrainZ = (int)round( x[2] / brainHz  );
                    
                    
                    Real tumourIC;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )
                    {
                        tumourIC  = 0.;
                    }
                    else
                    {
                        tumourIC  =  Tumor(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    
                    block(ix,iy,iz).phi = tumourIC * scale;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



void Glioma_RAT_UQ::_ic_rat_two_foci_tumor(Grid<W,B>& grid, int pID)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
    sprintf(dataFolder,"/home/baldesi/Glioma/RATGBM/source/Anatomy/F98/");
#endif
    
    sprintf(patientFolder, "%sM%02d/M%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s \n", patientFolder);
    
    sprintf(anatomy, "%s_gm.dat", patientFolder);
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%s_wm.dat", patientFolder);
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%s_csf.dat", patientFolder);
    MatrixD3D CSF(anatomy);
    sprintf(anatomy, "%s_mask.dat", patientFolder);
    MatrixD3D MASK(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.117;   // voxel spacing 117 µm, convert to mm -> L ~ 14 mm
    
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    std::cout<<"L="<<L<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(_DIM);
    vector<Real> tumor_ic2(_DIM);
    
     Real tumorRadius;
     Real tumorRadius2;

    _readInTumorPosition(tumor_ic);
    
   // IC are two identical sphere, which will fuse into elongated tumour
    if(pID == 1)
     {
       tumor_ic2[0] = 0.30;  
       tumor_ic2[1] = 0.54;
       tumor_ic2[2] = 0.48;
       tumorRadius  = 0.01;
       tumorRadius2 = 0.01;
     }
     
     // IC two sphere of different size, will shape two foci tumour
     if( pID ==2)
     {
       tumor_ic2[0] = 0.35;
       tumor_ic2[1] = 0.62;
       tumor_ic2[2] = 0.48;
       tumorRadius  = 0.011;
       tumorRadius2 = 0.007;
     }

    const Real smooth_sup  = 2;		// 2.suppor of smoothening, over how many gp to smooth
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const Real h = vInfo[0].h[0];  // make sure inference and propagation is run with same h for IC
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)round( x[0] / brainHx  );
                    int mappedBrainY = (int)round( x[1] / brainHy  );
                    int mappedBrainZ = (int)round( x[2] / brainHz  );
                    
                    
                    Real PGt, PWt, Pcsf, Pmask;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PGt   = 0.;
                        PWt   = 0.;
                        Pcsf  = 0.;
                        Pmask = 0.;
                    }
                    else
                    {
                        PGt  =  GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PWt  =  WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pcsf =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pmask = MASK(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    double all = PGt + PWt + Pcsf;
                    
                    if( all >  0.1 )
                    {
                        Pcsf = ( Pcsf > 0.1 ) ? 1. : Pcsf;  // threasholding to ensure hemisphere separations
                        block(ix,iy,iz).p_csf = Pcsf;
                        
                        if(Pcsf  < 1.)
                        {
                            if(PWt > 0.5)
                            {
                                block(ix,iy,iz).p_w    = 1.;//  / (PWt + PGt);
                                block(ix,iy,iz).p_g    = 0.;//  / (PWt + PGt);
                            }
                            else if (PGt > 0.5)
                            {
                                block(ix,iy,iz).p_w    = 0.;
                                block(ix,iy,iz).p_g    = 1.;
                            }
                        }
                        
                    }
                    
                    // fill the holes in the anatomy segmentations for the rats
                    all = block(ix,iy,iz).p_csf + block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    
                    if( (Pmask > 0.1 )&&( all< 0.1 )  )
                        block(ix,iy,iz).p_g = 1.;
                    
                    /* tumor foci 1*/
                    Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                    Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    Real psi = (dist - tumorRadius)*iw;
                    
                    if ((psi < -1)&& ((PGt>0.001) || (PWt >0.001)) )		// we are in tumor
                        block(ix,iy,iz).phi = 1.0;
                    else if(( (-1 <= psi) && (psi <= 1) )&& ((PGt>0) || (PWt >0)) )
                        block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                    else
                        block(ix,iy,iz).phi = 0.0;
                    
                    
                    /* tumor foci 2*/
                    p[0] = x[0] - tumor_ic2[0];
                    p[1] = x[1] - tumor_ic2[1];
                    p[2] = x[2] - tumor_ic2[2];
                    dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    psi = (dist - tumorRadius2)*iw;
                    
                    if ((psi < -1)&& ((PGt>0.001) || (PWt >0.001)) )		// we are in tumor
                        block(ix,iy,iz).phi += 1.0;
                    else if(( (-1 <= psi) && (psi <= 1) )&& ((PGt>0) || (PWt >0)) )
                        block(ix,iy,iz).phi += 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                    else
                        block(ix,iy,iz).phi += 0.0;
                    
                    /*scale tumour so max concentration is 1 to mimic tumour injection */
                    block(ix,iy,iz).phi = min(1.0, 4.*block(ix,iy,iz).phi ); //to ensure that max concentraiton is 1, for correct front propagation speed
                    block(ix,iy,iz).chi = Pmask;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}




void Glioma_RAT_UQ::_ic_rat_elongated_tumor(Grid<W,B>& grid, int pID)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
    sprintf(dataFolder,"/home/baldesi/Glioma/RATGBM/source/Anatomy/F98/");
#endif
    
    sprintf(patientFolder, "%sM%02d/M%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s \n", patientFolder);
    
    sprintf(anatomy, "%s_gm.dat", patientFolder);
    MatrixD3D GM(anatomy);
    sprintf(anatomy, "%s_wm.dat", patientFolder);
    MatrixD3D WM(anatomy);
    sprintf(anatomy, "%s_csf.dat", patientFolder);
    MatrixD3D CSF(anatomy);
    sprintf(anatomy, "%s_mask.dat", patientFolder);
    MatrixD3D MASK(anatomy);
    
    int brainSizeX = (int) GM.getSizeX();
    int brainSizeY = (int) GM.getSizeY();
    int brainSizeZ = (int) GM.getSizeZ();
    
    int brainSizeMax = max(brainSizeX, max(brainSizeY,brainSizeZ));
    L    = brainSizeMax * 0.117;   // voxel spacing 117 µm, convert to mm -> L ~ 14 mm
    
    std::cout<<"brainSizeX="<<brainSizeX<<" brainSizeY="<<brainSizeY<<" brainSizeZ="<<brainSizeZ<<std::endl;
    std::cout<<"L="<<L<<std::endl;
    
    double brainHx = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHy = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    double brainHz = 1.0 / ((double)(brainSizeMax)); // should be w.r.t. longest dimension for correct aspect ratio
    
    /* Tumor Set UP */
    vector<Real> tumor_ic(_DIM);
    vector<Real> tumor_ic2(_DIM);
    _readInTumorPosition(tumor_ic);
    
    // Tumour initial injection for F98: 3µL = 3 mm^3 --> volume corresponding to sphere with radius r=0.8947 mm
    
    const Real tumorRadius = 0.01;//0.8947/L ; // map to [0,1]^3 space
    const Real smooth_sup  = 2;		// 2.suppor of smoothening, over how many gp to smooth
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid.getBlockCollection()[info.blockID];
        
        const Real h = vInfo[0].h[0];  // make sure inference and propagation is run with same h for IC
        const Real iw = 1./(smooth_sup * h);   // width of smoothening => now it is over two grid points
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    /* Anatomy */
                    int mappedBrainX = (int)round( x[0] / brainHx  );
                    int mappedBrainY = (int)round( x[1] / brainHy  );
                    int mappedBrainZ = (int)round( x[2] / brainHz  );
                    
                    
                    Real PGt, PWt, Pcsf, Pmask;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PGt   = 0.;
                        PWt   = 0.;
                        Pcsf  = 0.;
                        Pmask = 0.;
                    }
                    else
                    {
                        PGt  =  GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PWt  =  WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pcsf =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pmask = MASK(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    double all = PGt + PWt + Pcsf;
                    
                    if( all >  0.1 )
                    {
                        Pcsf = ( Pcsf > 0.1 ) ? 1. : Pcsf;  // threasholding to ensure hemisphere separations
                        block(ix,iy,iz).p_csf = Pcsf;
                        
                        if(Pcsf  < 1.)
                        {
                            if(PWt > 0.5)
                            {
                                block(ix,iy,iz).p_w    = 1.;//  / (PWt + PGt);
                                block(ix,iy,iz).p_g    = 0.;//  / (PWt + PGt);
                            }
                            else if (PGt > 0.5)
                            {
                                block(ix,iy,iz).p_w    = 0.;
                                block(ix,iy,iz).p_g    = 1.;
                            }
                        }
                        
                    }
                    
                    // fill the holes in the anatomy segmentations for the rats
                    all = block(ix,iy,iz).p_csf + block(ix,iy,iz).p_w + block(ix,iy,iz).p_g;
                    
                    if( (Pmask > 0.1 )&&( all< 0.1 )  )
                        block(ix,iy,iz).p_g = 1.;
                    
                    /* tumor part 1:*/
                    Real p[3] = {x[0] - tumor_ic[0], x[1] - tumor_ic[1], x[2] - tumor_ic[2]};
                    Real dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                    Real psi = (dist - tumorRadius)*iw;
                    
                    if ((psi < -1)&& ((PGt>0.001) || (PWt >0.001)) )		// we are in tumor
                        block(ix,iy,iz).phi = 1.0;
                    else if(( (-1 <= psi) && (psi <= 1) )&& ((PGt>0) || (PWt >0)) )
                        block(ix,iy,iz).phi = 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                    else
                        block(ix,iy,iz).phi = 0.0;
                    
                    /* tumor part 2,3,4 */
                    for (int i=1; i<5; i++)
                    {
                        tumor_ic2[0] = tumor_ic[0] - i*0.01;
                        tumor_ic2[1] = tumor_ic[1] - i*0.01;
                        tumor_ic2[2] = tumor_ic[2] ;

                        p[0] = x[0] - tumor_ic2[0];
                        p[1] = x[1] - tumor_ic2[1];
                        p[2] = x[2] - tumor_ic2[2];

                        dist = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);    // distance of curent voxel from tumor center
                        psi = (dist - tumorRadius)*iw;
                        
                        if ((psi < -1)&& ((PGt>0.001) || (PWt >0.001)) )		// we are in tumor
                            block(ix,iy,iz).phi += 1.0;
                        else if(( (-1 <= psi) && (psi <= 1) )&& ((PGt>0) || (PWt >0)) )
                            block(ix,iy,iz).phi += 1.0 * 0.5 * (1 - psi - sin(M_PI*psi)/(M_PI));
                        else
                            block(ix,iy,iz).phi += 0.0;
                       
                    } 
                    
 
                    // scale so max tumor concentraiton is 1 to mimic tumour injection
                    block(ix,iy,iz).phi = min(1.0, 4.*block(ix,iy,iz).phi ); //to ensure that max concentraiton is 1, for correct front propagation speed
                    block(ix,iy,iz).chi = Pmask;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}




#pragma mark ReactionDiffusion
void Glioma_RAT_UQ::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    UpdateTumor                     <_DIM>  updateTumor(dt);
    
   blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}
    
    
void Glioma_RAT_UQ::_reactionDiffusionNecrosisStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, const double dt, const Real gamma)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_ReactionDiffusionNecrosisOperator<_DIM>  rhs(Dw,Dg,rho,gamma);
    UpdateTumorNecoris                      <_DIM>  updateTumor(dt);
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}


#pragma mark DumpingOutput
void Glioma_RAT_UQ:: _dump(int counter)
{
    
    if (parser("-vtk").asBool())
    {
        char filename[256];
        sprintf(filename,"%dD_%d_Rat_data%04d",_DIM, pID, counter);
        
        if( _DIM == 2)
        {
            IO_VTKNative<W,B, 2,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
        else
        {
            IO_VTKNative3D<W,B, 7,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
    }
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_RAT_UQ::_dumpUQoutput(double t)
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    double eps = hf*0.5;
    int time_point = (int) t;

    if(bVerbose) printf("bpd=%i, bs=%i, hf=%f,\n",blocksPerDimension,blockSize,hf);
    
    MatrixD3D tumor(gpd,gpd,gpd);
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
#pragma omp parallel for
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        double h = info.h[0];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    //mapped coordinates
                    int mx = (int)round( (x[0]) / hf  );
                    int my = (int)round( (x[1]) / hf  );
                    int mz = (int)round( (x[2]) / hf  );
                     
                    if( h < hf + eps)
                        tumor(mx,my,mz) = block(ix,iy,iz).phi;
		   else if(h < 2.*hf + eps)
                    {
                        for(int cz=0; cz<2; cz++)
                            for(int cy=0; cy<2; cy++)
                                for(int cx=0; cx<2; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    }
                    else if (h < 3.*hf + eps)
                    {
                        for(int cz=0; cz<3; cz++)
                            for(int cy=0; cy<3; cy++)
                                for(int cx=0; cx<3; cx++)
                                            tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                        
                    }
                    else
                    {
                        for(int cz=0; cz<4; cz++)
                            for(int cy=0; cy<4; cy++)
                                for(int cx=0; cx<4; cx++)
                                    tumor(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).phi;
                    			
			}
                }
        
        
    }
    
    char filename[256];
    sprintf(filename,"M_UQ_J%02d.dat",time_point);
    tumor.dump(filename);
    
}


void Glioma_RAT_UQ::run()
{
    
    bool bProfiler = 0;
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    /* scale to characteristic units */
    Real Dscale = 1./parser("-Dscale").asDouble();
    Dw = Dw/(L*L);  // rescale w.r.t. to characeteristic length
    Dg = Dscale*Dw;
    

    /* simulation set up */
    double tend         = parser("-Tend").asDouble(11.);
    double t            = 0;
    int iCounter        = 1;
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = 0.95 * h*h / ( 2.* _DIM * max(Dw, Dg) );
    if(bVerbose)  printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f, scale=%e \n", Dg, Dw, dt, rho,h, scale);
    
    /* set times of MRI scans */
    vector<int> timeOfscan ;
    timeOfscan.push_back(9);
    timeOfscan.push_back(11);
    timeOfscan.push_back(14);
    timeOfscan.push_back(16);
    std::vector<int>::iterator it = timeOfscan.begin();
    whenToWrite = *it;
    
    // since inpute file contains tumour at day 9, save just day 11
    if(ICtype == 1)
        t = 9;
    
    
    
    /* initial refinement & compression */
    if( (whenToRefine > 0.) && (bAdaptivity) )
    {
        Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
        Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
    }
    
    while (t <= tend)
    {
        _reactionDiffusionStep(boundaryInfo, nParallelGranularity, Dw, Dg, rho, dt);
        t                   += dt   ;
        numberOfIterations  ++      ;
        
        
        // refinment & compression
        if (( t >= ((double)(whenToRefine)) ) && (bAdaptivity) )
        {
            Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
            Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
           
            whenToRefine = whenToRefine+whenToRefineOffset;
        }
        
        if ( t >= ((double)(whenToWrite)) )
        {
            _dumpUQoutput(whenToWrite);
            _dump((int) whenToWrite );
            
            ++it;
            whenToWrite = *it;
            if(bVerbose) printf("Dumping data at time t=%f\n", t);
            
        }
    }

  // save the last point if omited because of too big time step
  if(t > whenToWrite)
  {
            if(bVerbose) printf("Dumping data at time t=%f\n", t);
      _dumpUQoutput(whenToWrite);
      _dump((int) whenToWrite );
  }    


    // Refine final state & dump for UQ Likelihood
//    if(bAdaptivity)
//        Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);

    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
