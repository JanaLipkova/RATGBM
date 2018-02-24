//
//  Glioma_RAT_preprocessing.h
//  RATGBM_xcode
//
//  Created by Lipkova on 14/02/18.
//  Copyright (c) 2018 Lipkova. All rights reserved.
//

// 1) Read in anatomy, mask and tumour segmentations
// 2) compute center of mass of the tumour segmentation + volume -> save to file
// 3) save tumour segm + mask into binary file

#pragma once
#include "Glioma_Types.h"


class Glioma_RAT_preprocessing: public Glioma
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
    bool                                    isDone;
    bool                                    bAdaptivity;
    bool                                    bVerbose;
    int                                     pID;
    Real                                    VT1,VT2;
    
    static void _ic(Grid<W,B>& grid, int pID);
    void        _computeTumourProperties(bool bDumbIC2file);
    void        _computeEnclosingSphere(Grid<W,B>& grid);
    void        _readInTumorPosition(vector<Real>& tumorIC );
    void		_dump(int counter);
    void        _dump2binary(int day);
    void        _readInTumourSegmentation(Grid<W,B>& grid, int pID, int day);

    
    
public:
    Glioma_RAT_preprocessing(int argc, const char ** argv);
    ~Glioma_RAT_preprocessing();
    void run();
    
};

