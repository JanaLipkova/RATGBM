set EXTRA_INC=$HOME/Glioma/sourcestochasticF
g++ -O3 -w -I$EXTRA_INC ./MRAGcore/*.cpp
g++ -O3 -w -I$EXTRA_INC ./RD/SRD.cpp


mpi-inc = ${MPI_ROOT}/include/
mpi-lib = ${MPI_ROOT}/lib/
		
TBB_INC_DIR=/cluster/home/mavt/lipkovaj/LIB/tbb40_20120613oss/include
TBB_LIB_DIR=/cluster/home/mavt/lipkovaj/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release
	
VTK_INC_DIR=/cluster/home/mavt/lipkovaj/LIB/myVTK/include/vtk-5.2
VTK_LIB_DIR=/cluster/home/mavt/lipkovaj/LIB/myVTK/lib/vtk-5.2
	
export LD_LIBRARY_PATH=/cluster/home/mavt/lipkovaj/LIB/myVTK/lib/vtk-5.2:$LD_LIBRARY_PATH
	
CPPFLAGS+= -I$(TBB_INC_DIR) -I$(VTK_INC_DIR)
	
LIBS += \
	-L$(TBB_LIB_DIR) \
	-ltbb \
	-ltbbmalloc \
	-L$(VTK_LIB_DIR) \
	-lvtkViews \
	-lvtkInfovis \
	-lvtkWidgets \
	-lvtkHybrid \
	-lvtkRendering \
	-lvtkGraphics \
	-lvtkverdict \
	-lvtkImaging \
	-lvtkftgl \
	-lvtkfreetype \
	-lvtkIO \
	-lvtkFiltering \
	-lvtkCommon \
	-lm \
	-lvtkDICOMParser \
	-lvtkmetaio \
	-lvtksqlite \
	-lvtkpng \
	-lvtktiff \
	-lvtkjpeg \
	-lvtkexpat \
	-lvtksys \
	-lvtkexoIIc \
	-lvtkNetCDF \
	-lvtklibxml2 \
	-lvtkzlib \
	-lpthread \
	-ldl \
	$(OPENMP_FLAG)
	
	
##################	

set EXTRA_INC=$HOME/Glioma/sourcestochasticF
g++ -O3 -w -I$EXTRA_INC ./MRAGcore/*.cpp
g++ -O3 -w -I$EXTRA_INC ./RD/SRD.cpp

g++ SRD.o MRAGBoundaryBlockInfo.o MRAGProfiler.o MRAGWavelets_StaticData.o -o brain $(LIBS)




