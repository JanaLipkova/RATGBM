//
//  Glioma_RAT_IC_Approximation.h
//  RATGBM_xcode
//
//  Created by Lipkova on 24/02/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//


#pragma once
#include "Glioma_Types.h"


class Glioma_RAT_IC_Approximation: public Glioma
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
    
    
    static void _ic(Grid<W,B>& grid, int pID);
    static void _readInTumourSegmentation(Grid<W,B>& grid, int pID, int day);
    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt);
    void		_dump(int counter);
    void        _dump2binary(int counter);
    void        _normaliseTumour();

public:
    Glioma_RAT_IC_Approximation(int argc, const char ** argv);
    ~Glioma_RAT_IC_Approximation();
    void run();
    
};

