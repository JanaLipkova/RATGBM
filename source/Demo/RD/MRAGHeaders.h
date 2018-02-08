#include "../../MRAGcore/MRAGCommon.h"
#include "../../MRAGcore/MRAGEnvironment.h"

#define _WAVELET_TYPE Wavelets_AverageInterp3rdOrder
#ifdef _MRAG_TBB
    static const int nThreads = _MRAG_TBB_NTHREADS_HINT ;
#else
    static const int nThreads = 1 ;
#endif
static const int blockSize = 16 ;
static const int blocksPerDimension = 128/32;
static const int maxLevel = 3;



double t = 0;
double dt = 0;
static const int maxStencil[2][3] = {
-1, -1, -1,
+2, +2, +2
};



#if _MRAG_OS == _MRAG_OS_APPLE
#ifdef _MRAG_GLUT_VIZ
#include "GLUT/glut.h"
#endif
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#define _USE_MATH_DEFINES
#ifdef _MRAG_GLUT_VIZ
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
#ifdef _MRAG_TBB
    #include "../../MRAGmultithreading/MRAGBlockProcessing_TBB.h"
#endif

#include "../../MRAGio/MRAG_IO_Native.h"
#include "../../MRAGio/MRAG_IO_VTK.h"
