/*
 *  CompressibleFlow.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/8/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#include "MRAGcore/MRAGEnvironment.h"

//#define __DEVICE_EMULATION__ 1
#include <cutil.h>
#include <cuda_runtime_api.h>

#if _MRAG_OS == _MRAG_OS_APPLE
#include "GLUT/glut.h"
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#include <stdlib.h>
#include <stdio.h>
#include "GL/glew.h"
#include "GL/glut.h"
#define _USE_MATH_DEFINES
#endif

#include <limits>

using namespace std;

#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGBlockCollection.h"
#include "MRAGcore/MRAGRefiner.h"
#include "MRAGcore/MRAGProfiler.h"
#include "MRAGcore/MRAGBlockLab.h"

#include "../MRAGmultithreading/MRAGBlockProcessing_CUDA.h"
#include "../MRAGmultithreading/MRAGBlockCollection_CUDA.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_TBB.h"

#include "../MRAGvisual/GridViewer.h"
#include "../MRAGvisual/MRAGVisualTypes.h"

#include "../MRAGscience/MRAGScienceCore.h"
#include "../MRAGscience/MRAGSpaceTimeSorter.h"
#include "../MRAGio/MRAG_IO_Native.h"

//typedef MRAG::Wavelets_Interp2ndOrder W;
typedef MRAG::Wavelets_Interp4thOrder W;


#include "CompressibleFlow.h"
#include "CompressibleFlowTypes.h"

#ifdef _MRAG_TBB
typedef Multithreading::BlockProcessing_TBB<B> BlockProcessingCPU;
#else
typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessingCPU;
#endif

#include "Informal_DiffusionTest_CUDA.h"

using namespace MRAG;


GridViewer viewer(false, true);
DiffusionTestCUDA simulation; 


RGBA convertToRGBA(FluidElement& p)
{
	RGBA c(max(0.0, 2.0*(p.phi+2.223)/4.445 - 0.5), 1.0 - 2.0*fabs((p.phi+2.223)/4.445 - 0.5),max(0.0, 1.0-2.0*(p.phi+2.223)/4.445),0);
	//RGBA c(max(0.0, 2.0*(p.rho+2.223)/4.445 - 0.5),  1.0 - 2.0*fabs((p.rho+2.223)/4.445 - 0.5),min(1.0, 0.0+2.0*(p.rho+2.223)/4.445),0);
	return c;
}

CompressibleFlow::CompressibleFlow(int argc, char ** argv): 
	time(0)
{
	CUT_DEVICE_INIT(argc, argv);
	
	simulation.setup();
}

void CompressibleFlow::Step()
{
	simulation.step();
}

void CompressibleFlow::Render()
{
	
	viewer.drawContent<W,B>(*simulation.grid, simulation.grid->getBlockCollection());
	viewer.drawSketch(*simulation.grid, false);

}

double CompressibleFlow::getTime()
{
	return time;
}

