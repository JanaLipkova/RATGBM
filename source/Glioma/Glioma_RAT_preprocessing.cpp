//
//  Glioma_RAT_preprocessing.cpp
//  RATGBM_xcode
//
//  Created by Lipkova on 14/02/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#include "Glioma_RAT_preprocessing.h"

static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_RAT_preprocessing::Glioma_RAT_preprocessing(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////            RAT GLIOMA PREPROCESSING          ////////////////\n");
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n", nThreads, blockSize, "w", blocksPerDimension, maxLevel);
    
    refiner		= new Refiner_SpaceExtension(resJump,maxLevel);
    compressor	= new Compressor(resJump);
    Environment::setup();
    
    grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
    grid->setCompressor(compressor);
    grid->setRefiner(refiner);
    stSorter.connect(*grid);
    
    pID =  parser("-pID").asInt();
    _ic(*grid, pID);
    
    isDone              = false;
}

Glioma_RAT_preprocessing::~Glioma_RAT_preprocessing()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
/* Initisalise tumour as a point sources, i.e. small smooth sphere
 1) read in anatomies - rescaled to [0,1]^3
 2) read in tumor center of mass + initialize tumor around
 3) set length of brain */
void Glioma_RAT_preprocessing::_ic(Grid<W,B>& grid, int pID)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
    sprintf(dataFolder,"../../Anatmoy/");
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
    sprintf(anatomy, "%sJ09-T2w-iso-crop.dat", patientFolder);
    MatrixD3D T2w(anatomy);
    sprintf(anatomy, "%sJ09-DCE-iso-crop.dat", patientFolder);
    MatrixD3D T1w(anatomy);
    
    
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

  // for synthethic data
    double x1, x2;
    const float alpha = 0.02;

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
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    
                    Real PGt, PWt, Pcsf, PT2w, PT1w, Pmask;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PGt   = 0.;
                        PWt   = 0.;
                        Pcsf  = 0.;
                        PT2w  = 0.;
                        PT1w  = 0.;
                        Pmask = 0.;
                    }
                    else
                    {
                        PGt   =  GM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PWt   =  WM(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pcsf  =  CSF(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PT2w  =  T2w(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PT1w  =  T1w(mappedBrainX,mappedBrainY,mappedBrainZ);
                        Pmask =  MASK(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    double all = PGt + PWt + Pcsf;
                    
                    if( (all > 0.) || (Pmask > 0.) )
                    {
                        Pcsf = ( Pcsf > 0.1 ) ? 1. : Pcsf;  // threasholding to ensure hemisphere separations
                        block(ix,iy,iz).p_csf = Pcsf;
                        
                        if(Pcsf  < 1.)
                        {
                            block(ix,iy,iz).p_w    = PWt  / (PWt + PGt);
                            block(ix,iy,iz).p_g    = PGt  / (PWt + PGt);
                        }
                    }
                    
                    if(pID > 0)
                    {   
                        block(ix,iy,iz).t1bc = (PT1w > 0.01) ? 1. : 0. ;  // T1w tumour segmentations
                        block(ix,iy,iz).t2bc = (PT2w > 0.01) ? 1. : 0.;   // T2w tumour segmentations
                    }
                    else
                    {
                       
                       if (PT1w > 0.01)
                       {
                          /*
                          if(addedNoise)
                          {
                          // Box-Muller
			   double u1 = drand48();
                           double u2 = drand48();

                           x1 = sqrt( - 2.0*log( u1 ) ) * cos( 2.0*M_PI*u2 );
                           x2 = sqrt( - 2.0*log( u1 ) ) * sin( 2.0*M_PI*u2 );
                           
			   double ucT1 = max( 0., 0.7  + alpha * x1);
                           double ucT2 = max( 0., 0.25 + alpha * x2);
                           }*/

                           double ucT1 = 0.7;
                           double ucT2 = 0.25;
                           block(ix,iy,iz).t1bc = (PT1w >= ucT1) ? 1. : 0. ;  // T1w tumour segmentations
                           block(ix,iy,iz).t2bc = (PT2w >= ucT2) ? 1. : 0. ;  // T2w tumour segmentations
                       }
			
		   }
                    
                    block(ix,iy,iz).chi = Pmask;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}



#pragma mark Diagnostic
void Glioma_RAT_preprocessing::_computeTumourProperties(bool bDumbIC2file)
{
    Real massT1 = 0.;
    Real massT2 = 0.;
    vector<Real> cmT1(3,0.);
    vector<Real> cmT2(3,0.);

    double h;
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        h = info.h[0];

        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    double x[3];
                    info.pos(x, ix, iy, iz);
                    
                    massT1  += block(ix,iy,iz).t1bc;
                    cmT1[0] += x[0] * block(ix,iy,iz).t1bc;
                    cmT1[1] += x[1] * block(ix,iy,iz).t1bc;
                    cmT1[2] += x[2] * block(ix,iy,iz).t1bc;

                    massT2  += block(ix,iy,iz).t2bc;
                    cmT2[0] += x[0] * block(ix,iy,iz).t2bc;
                    cmT2[1] += x[1] * block(ix,iy,iz).t2bc;
                    cmT2[2] += x[2] * block(ix,iy,iz).t2bc;
                }
        
    }
    
    cmT1[0] = cmT1[0] / massT1;
    cmT1[1] = cmT1[1] / massT1;
    cmT1[2] = cmT1[2] / massT1;
    
    
    cmT2[0] = cmT2[0] / massT2;
    cmT2[1] = cmT2[1] / massT2;
    cmT2[2] = cmT2[2] / massT2;

    
    printf("Center of mass T1: cx=%f, cy=%f, cz=%f \n", cmT1[0], cmT1[1], cmT1[2]);
    printf("Center of mass T2: cx=%f, cy=%f, cz=%f \n", cmT2[0], cmT2[1], cmT2[2]);

    Real h3 = h*h*h;
    Real VT1 = massT1 * h3;
    Real VT2 = massT2 * h3;
    
    printf("VT1 = %f, VT2 = %f \n", VT1 , VT2);
    
    if(bDumbIC2file)
    {
        Real r3 =VT1 * 3. / ( 4. * M_PI);
        Real r = pow(r3,1./3.);
        printf("Prior range: \n");
        printf("icx = [ %f, %f] \n", cmT1[0] - r, cmT1[0] + r);
        printf("icy = [ %f, %f] \n", cmT1[1] - r, cmT1[1] + r);
        printf("icz = [ %f, %f] \n", cmT1[2] - r, cmT1[2] + r);
        
        FILE * pFile;
        float buffer[3] = { cmT1[0] , cmT1[1] , cmT1[2] };
        pFile = fopen ("HGG_TumorIC.bin", "wb");
        fwrite (buffer , sizeof(float), sizeof(buffer), pFile);
        fclose (pFile);
    }
}


void Glioma_RAT_preprocessing:: _readInTumourSegmentation(Grid<W,B>& grid, int pID, int day)
{
    char dataFolder   [200];
    char patientFolder[200];
    char anatomy      [200];
    
#ifdef LRZ_CLUSTER
    sprintf(dataFolder,"/home/hpc/txh01/di49zin/GliomaAdvance/RATGBM/source/Anatomy/F98/");
#else
    sprintf(dataFolder,"../../Anatmoy/");
#endif
    
    sprintf(patientFolder, "%sM%02d/M%02d",dataFolder,pID,pID);
    printf("Reading anatomy from: %s \n", patientFolder);
    
    sprintf(anatomy, "%sJ%02d-T2w-iso-crop.dat", patientFolder,day);
    MatrixD3D T2w(anatomy);
    sprintf(anatomy, "%sJ%02d-DCE-iso-crop.dat", patientFolder,day);
    MatrixD3D T1w(anatomy);
    
    
    int brainSizeX = (int) T2w.getSizeX();
    int brainSizeY = (int) T2w.getSizeY();
    int brainSizeZ = (int) T2w.getSizeZ();
    
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
                    int mappedBrainX = (int)floor( x[0] / brainHx  );
                    int mappedBrainY = (int)floor( x[1] / brainHy  );
                    int mappedBrainZ = (int)floor( x[2] / brainHz  );
                    
                    
                    Real PT2w, PT1w;
                    
                    if ( (mappedBrainX < 0 || mappedBrainX >= brainSizeX) || (mappedBrainY < 0 || mappedBrainY >= brainSizeY) || (mappedBrainZ < 0 || mappedBrainZ >= brainSizeZ) )                    {
                        PT2w  = 0.;
                        PT1w  = 0.;
                    }
                    else
                    {
                        PT2w  =  T2w(mappedBrainX,mappedBrainY,mappedBrainZ);
                        PT1w  =  T1w(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    if(pID > 0)
                    {
                        block(ix,iy,iz).t1bc = (PT1w > 0.01) ? 1. : 0. ;  // T1w tumour segmentations
                        block(ix,iy,iz).t2bc = (PT2w > 0.01) ? 1. : 0.;   // T2w tumour segmentations
                    }
                    else
                    {
                        Real ucT1 = 0.7;
                        Real ucT2 = 0.25;
                        
                        block(ix,iy,iz).t1bc = (PT1w >= ucT1) ? 1. : 0. ;  // T1w tumour segmentations
                        block(ix,iy,iz).t2bc = (PT2w >= ucT2) ? 1. : 0. ;  // T2w tumour segmentations
                    }
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }

    
    
}

#pragma mark DumpingOutput
void Glioma_RAT_preprocessing:: _dump(int counter)
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
            IO_VTKNative3D<W,B, 9,0 > vtkdumper2;
            vtkdumper2.Write(*grid, grid->getBoundaryInfo(), filename);
        }
    }
}

/* Dump output for UQ likelihood. Requirements:
 - dump at the uniform finest resolution
 - use 3D Matrix structure to dump data in binary format
 - assume 3D simulation */
void Glioma_RAT_preprocessing::_dump2binary(int day)
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    double eps = hf*0.1;
    
    if(bVerbose) printf("bpd=%i, bs=%i, hf=%f,\n",blocksPerDimension,blockSize,hf);
    
    MatrixD3D tumorT1(gpd,gpd,gpd);
    MatrixD3D tumorT2(gpd,gpd,gpd);
    MatrixD3D brainMask(gpd,gpd,gpd);


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
                    int mx = (int)floor( (x[0]) / hf  );
                    int my = (int)floor( (x[1]) / hf  );
                    int mz = (int)floor( (x[2]) / hf  );
                    
                    if(h < 2.*hf - eps)
                    {
                        tumorT1(mx,my,mz)   = block(ix,iy,iz).t1bc;
                        tumorT2(mx,my,mz)   = block(ix,iy,iz).t2bc;
                        brainMask(mx,my,mz) = block(ix,iy,iz).chi;

                    }
                    else if(h < 3.*hf - eps)
                    {
                        for(int cz=0; cz<2; cz++)
                            for(int cy=0; cy<2; cy++)
                                for(int cx=0; cx<2; cx++)
                                {
                                    tumorT1(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).t1bc;
                                    tumorT2(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).t2bc;
                                    brainMask(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).chi;
                                }
                        
                    }
                    else if (h < 4.*hf - eps)
                    {
                        for(int cz=0; cz<3; cz++)
                            for(int cy=0; cy<3; cy++)
                                for(int cx=0; cx<3; cx++)
                                {
                                    tumorT1(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).t1bc;
                                    tumorT2(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).t2bc;
                                    brainMask(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).chi;
                                }
                        
                    }
                    else
                    {
                        for(int cz=0; cz<4; cz++)
                            for(int cy=0; cy<4; cy++)
                                for(int cx=0; cx<4; cx++)
                                {
                                    tumorT1(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).t1bc;
                                    tumorT2(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).t2bc;
                                    brainMask(mx+cx,my+cy,mz+cz) = block(ix,iy,iz).chi;

                                }
                    }
                }
        
    }
    
    char filename[256];
    
    sprintf(filename,"T1w_D%02d.dat",day);
    tumorT1.dump(filename);
    
    sprintf(filename,"T2w_D%02d.dat",day);
    tumorT2.dump(filename);
    
    sprintf(filename,"Mask.dat");
    brainMask.dump(filename);
    
}


void Glioma_RAT_preprocessing::run()
{
    
    printf("run start \n");
    
    vector<int> timeOfscan ;
    timeOfscan.push_back(9);
    timeOfscan.push_back(11);
    timeOfscan.push_back(14);
    timeOfscan.push_back(16);
    bool bDumbIC2file;
    
    for (std::vector<int>::iterator it = timeOfscan.begin(); it  !=timeOfscan.end() ; ++it)
    {
        std::cout<<"Tumour statistics at the day="<< *it << std::endl;
        
        if(*it == 9)
        {
            bDumbIC2file = 1;
            _computeTumourProperties(bDumbIC2file);  // compute center of mass and tumour volume
            _dump(*it);
            _dump2binary(*it);
        }
        else
        {
            bDumbIC2file = 0;
            _readInTumourSegmentation(*grid, pID, *it);
            _computeTumourProperties(bDumbIC2file);  // compute center of mass and tumour volume
            _dump(*it);
            _dump2binary(*it);
        }
    }
        
    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
