//
//  Glioma_HG_UQ.h
//  GliomaXcode
//
//  Created by Lipkova on 10/06/15.
//  Copyright (c) 2015 Lipkova. All rights reserved.
//

#ifndef Glioma_HG_UQ__
#define Glioma_HG_UQ__

#pragma once
#include "Glioma_Types.h"


class Glioma_HG_UQ: public Glioma
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
    
    static void _ic_SubjectBrain(Grid<W,B>& grid);
    static void _ic_AtlasBrain(Grid<W,B>& grid);
    static void _ic_PatientCase(Grid<W,B>& grid, int pID);

    static void	_readInBrainWebAnatomy(vector<float>& tissue, FILE* fp, int DataSize, int threshold );
    static void _selectBrainWebAnatomy(vector<float>& GreyTissueData, vector<float>& WhiteTissueData, vector<float> & CsfData, int dataSize);
    static void _readInTumorPosition(vector<Real>& tumorIC);

    void        _reactionDiffusionStep(BoundaryInfo* boundaryInfo, const int nParallelGranularity, const Real Dw, const Real Dg, const Real rho, double dt);
    void		_dump(int counter);
    void        _dumpUQoutput(Grid<W,B>& grid);
    void        _dumpBrainPoints(Grid<W,B>& grid );
    void        _dumpSubBrainPoints(Grid<W,B>& grid );


public:
    Glioma_HG_UQ(int argc, const char ** argv);
    ~Glioma_HG_UQ();
    void run();
    
};

#endif /* defined(Glioma_HG_UQ__) */
