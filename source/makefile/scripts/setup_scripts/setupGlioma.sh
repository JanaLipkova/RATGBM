module unload gcc
module load gcc/4.4.6

export LD_LIBRARY_PATH=/cluster/work/infk/wvanrees/apps/TBB/tbb41_20120718oss/build/linux_intel64_gcc_cc4.7.0_libc2.12_kernel2.6.32_release/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/cluster/home/mavt/lipkovaj/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/work/infk/cconti/hypre-2.10.0b/src/hypre/lib/:$LD_LIBRARY_PATH

