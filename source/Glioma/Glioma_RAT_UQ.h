//
//  Glioma_RAT_UQ.h
//  
//
//  Created by Lipkova on 08/02/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//
//

#pragma once
#include "Glioma_Types.h"


class Glioma_RAT_UQ: public Glioma
{
private:
    Grid<W, B>								* grid;
    BlockProcessing							blockProcessing;
    Refiner_SpaceExtension					*refiner;
    Compressor								*compressor;
    BlockFWT<W, B, RD_Projector_Wavelets>	blockfwt;			// refinment based on single channel
    SpaceTimeSorter							stSorter;
    Profiler								profiler;
    ArgumentParser							parser;
    IO_VTK< W, B, RD_Projector_VTK >		vtk;
    BlockLab< B >							lab;
    int										numberOfIterations;
    double                                  whenToWrite;
    double                                  whenToWriteOffset;
    bool                                    isDone;
    bool                                    bAdaptivity;
    bool                                    bVerbose;
    int                                     pID;
    int                                     ICtype;
    int                                     ModelType;
    double                                  scale;
    
    
    static void _ic_rat_point_tumor(*grid, pID);
    static void _ic_rat_tumour(*grid, pID);
    static void _readInTumorPosition(vector<Real>& tumorIC);
    void        _rescale_init_tumour(scale);
    
    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt);
    void		_dump(int counter);
    void        _dumpUQoutput(Grid<W,B>& grid);
    
    
public:
    Glioma_RAT_UQ(int argc, const char ** argv);
    ~Glioma_RAT_UQ();
    void run();
    
};

#endif /* defined(Glioma_HG_UQ__) */

