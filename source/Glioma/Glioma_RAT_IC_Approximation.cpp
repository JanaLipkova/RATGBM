//
//  Glioma_RAT_IC_Approximation.cpp
//  RATGBM_xcode
//
//  Created by Lipkova on 24/02/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

#include "Glioma_RAT_IC_Approximation.h"


static int maxStencil[2][3] = {
    -1, -1, -1,
    +2, +2, +2
};

Glioma_RAT_IC_Approximation::Glioma_RAT_IC_Approximation(int argc, const char ** argv): parser(argc, argv)
{
    bVerbose = parser("-verbose").asBool();
    
    if(bVerbose) printf("////////////////////////////////////////////////////////////////////////////////\n");
    if(bVerbose) printf("//////////////////         RAT GLIOMA IC APPROXIMATION          ////////////////\n");
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
    _readInTumourSegmentation(*grid, pID, 9);
    
    isDone              = false;
    
    whenToWriteOffset	= parser("-dumpfreq").asDouble();
    whenToWrite			= parser("-dumpstart").asDouble();
    numberOfIterations	= 0;
}

Glioma_RAT_IC_Approximation::~Glioma_RAT_IC_Approximation()
{
    std::cout << "------Adios muchachos------" << std::endl;
}


#pragma mark InitialConditions
/* Initisalise tumour as a point sources, i.e. small smooth sphere
 1) read in anatomies - rescaled to [0,1]^3
 2) read in tumor center of mass + initialize tumor around
 3) set length of brain */
void Glioma_RAT_IC_Approximation::_ic(Grid<W,B>& grid, int pID)
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
                        Pmask =  MASK(mappedBrainX,mappedBrainY,mappedBrainZ);
                    }
                    
                    double all = PGt + PWt + Pcsf;
                    
                    if( all > 0.1 )
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


void Glioma_RAT_IC_Approximation:: _readInTumourSegmentation(Grid<W,B>& grid, int pID, int day)
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
     
    const Real ucT1 = (pID == 1) ? 0.75 : 0.7;
    const Real ucT2 = 0.25;
  
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

                  if( pID > 10)
                  {                    
                    block(ix,iy,iz).t1bc = (PT1w > 0.01) ? 1. : 0. ;  // T1w tumour segmentations
                    block(ix,iy,iz).t2bc = (PT2w > 0.01) ? 1. : 0.;   // T2w tumour segmentations
                  }
                  else
		  { 	
		      block(ix,iy,iz).t1bc = (PT1w > ucT1)  ? 1. : 0.; // Threshodl solution
	              block(ix,iy,iz).t2bc = (PT2w > ucT2) ? 1. : 0.;

	          }
         
                    
                    if( (block(ix,iy,iz).t1bc > 0.5) && ( block(ix,iy,iz).t2bc == 0.) )
                        block(ix,iy,iz).phi = 1;
                    else
                        block(ix,iy,iz).phi = 0.5 * block(ix,iy,iz).t1bc + 0.5 * block(ix,iy,iz).t2bc;
                    
                }
        
        grid.getBlockCollection().release(info.blockID);
        
    }
}




#pragma mark ReactionDiffusion
void Glioma_RAT_IC_Approximation::_reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt)
{
    
    vector<BlockInfo> vInfo				= grid->getBlocksInfo();
    const BlockCollection<B>& collecton = grid->getBlockCollection();
    
    Glioma_ReactionDiffusionOperator<_DIM>  rhs(Dw,Dg,rho);
    UpdateTumor                     <_DIM>  updateTumor(dt);
    
    blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, rhs);
    BlockProcessing::process(vInfo, collecton, updateTumor, nParallelGranularity);
}



void Glioma_RAT_IC_Approximation::_normaliseTumour()
{
    Real maxTum = 0.;
    
    vector<BlockInfo> vInfo = grid->getBlocksInfo();
    
    
    // get max concentration
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    maxTum = max(maxTum, block(ix,iy,iz).phi);
    }
    
    Real iscale = 1./maxTum;
    for(int i=0; i<vInfo.size(); i++)
    {
        BlockInfo& info = vInfo[i];
        B& block = grid->getBlockCollection()[info.blockID];
        
        for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                    block(ix,iy,iz).phi = block(ix,iy,iz).phi * iscale;
    }
    
}


#pragma mark DumpingOutput
void Glioma_RAT_IC_Approximation:: _dump(int counter)
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
void Glioma_RAT_IC_Approximation::_dump2binary()
{
    int gpd = blocksPerDimension * blockSize;
    double hf  = 1./gpd;
    double eps = hf*0.1;
    
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
                    {
                        tumor(mx,my,mz)   = block(ix,iy,iz).phi;
                    }
                    else if(h < 2.*hf + eps)
                    {
                        for(int cz=0; cz<2; cz++)
                            for(int cy=0; cy<2; cy++)
                                for(int cx=0; cx<2; cx++)
                                {
                                    tumor(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).phi;
                                }
                        
                    }
                    else if (h < 3.*hf + eps)
                    {
                        for(int cz=0; cz<3; cz++)
                            for(int cy=0; cy<3; cy++)
                                for(int cx=0; cx<3; cx++)
                                {
                                    tumor(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).phi;
                                }
                        
                    }
                    else
                    {
                        for(int cz=0; cz<4; cz++)
                            for(int cy=0; cy<4; cy++)
                                for(int cx=0; cx<4; cx++)
                                {
                                    tumor(mx+cx,my+cy,mz+cz)   = block(ix,iy,iz).phi;
                                }
                    }
                }
        
        
    }
    
    char filename[256];
    
    sprintf(filename,"M%02d_TumourIC.dat",pID);
    tumor.dump(filename);

}


void Glioma_RAT_IC_Approximation::run()
{
    const int nParallelGranularity	= (grid->getBlocksInfo().size()<=8 ? 1 : 4);
    BoundaryInfo* boundaryInfo		= &grid->getBoundaryInfo();
    
    /* read in case specific parameters*/
    ifstream mydata("HGG_InputParameters.txt");
    Real Dg, Dw, rho;
    double tend = parser("-Tend").asDouble(11.);
    int iCounter = 1;
    
    if (mydata.is_open())
    {
        mydata >> Dw;
        mydata >> rho;
        mydata.close();
    }

    Real Dscale = 1./parser("-Dscale").asDouble();
    Dw = Dw/(L*L);  // rescale w.r.t. to characeteristic length
    Dg = Dscale*Dw;
    
    rho = 0.;
    
    /* simulation set up */
    double t			= 0.0;
    double h            = 1./(blockSize*blocksPerDimension);
    double dt           = 0.95 * h*h / ( 2.* _DIM * max(Dw, Dg) );
    if(bVerbose)  printf("Dg=%e, Dw=%e, dt= %f, rho=%f , h=%f\n", Dg, Dw, dt, rho,h);
    
    /* initial refinement & compression */
    if( (whenToWrite > 0.) && (bAdaptivity) )
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
        if ( t >= (double)whenToWrite )
        {
            if(bAdaptivity)
            {
                Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
                Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
            }
            
            _dump(iCounter);
            iCounter++;
            whenToWrite = whenToWrite+whenToWriteOffset;
            if(bVerbose) printf("Dumping data at time t=%f\n", t);

        }

    }
    
    _normaliseTumour();
    _dump(iCounter);
    _dump2binary();
    
    if(bVerbose) profiler.printSummary();
    if(bVerbose) printf("**** Dumping done\n");
    if(bVerbose) printf("\n\n Run Finished \n\n");
}
