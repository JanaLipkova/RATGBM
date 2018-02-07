#include "../../MRAGcore/MRAGCommon.h"
#include "../../MRAGcore/MRAGEnvironment.h"

#define _WAVELET_TYPE Wavelets_AverageInterp5thOrder
static const int nThreads = 2 ;
static const int blockSize = 16;
static const int blocksPerDimension = 128/32;
static const int maxLevel=2;



double t = 0;
double dt = 0;
static const int maxStencil[2][3] = {
	-3, -3, -3,
	+4, +4, +4
};

const double refinement_tolerance = 1e-3;
const double compression_tolerance = 2e-4;



#if _MRAG_OS == _MRAG_OS_APPLE
#ifdef _MRAG_GLUT_VIZ
#include "GLUT/glut.h"
#endif
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#define _USE_MATH_DEFINES
#ifdef _MRAG_GLUT_VIZ
#include <stdlib.h>
#include "GL/glew.h"
#include "GL/glut.h"
#endif
#endif

#undef min
#undef max

#include "../../MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "../../MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "../../MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "../../MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "../../MRAGcore/MRAGWavelets_Haar.h"
#include "../../MRAGcore/MRAGrid.h"
#include "../../MRAGcore/MRAGRefiner.h"
#include "../../MRAGcore/MRAGCompressor.h"
#include "../../MRAGcore/MRAGBlockLab.h"
#include "../../MRAGcore/MRAGBlockFWT.h"

#ifdef _MRAG_GLUT_VIZ
#include "../../MRAGvisual/GridViewer.h"
#endif

#include "../../MRAGscience/MRAGScienceCore.h"
#include "../../MRAGscience/MRAGAutomaticRefiner.h"
#include "../../MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "../../MRAGscience/MRAGSpaceTimeSorter.h"
#include "../../MRAGscience/MRAGRefiner_SpaceExtension.h"

#include "../../MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "../../MRAGmultithreading/MRAGBlockProcessing_TBB.h"

#include "../../MRAGio/MRAG_IO_Native.h"
#include "../../MRAGio/MRAG_IO_VTK.h"
