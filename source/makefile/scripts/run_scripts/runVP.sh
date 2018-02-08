#===============================
#  Script for GLIOMA TMCMC
#  * run brain simulaiton
#===============================


export LD_LIBRARY_PATH=/cluster/home/mavt/lipkovaj/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/:$LD_LIBRARY_PATH

N=16
export OMP_NUM_THREADS=$N
module unload gcc
module load gcc/4.4.6
module load open_mpi
# brain simulation set up
program=brain
vtk=1
vp=0
dumpfreq=50
verbose=1
anatomy=VPvisualisation

bsub -W 00:30 -n $N ./$program -nthreads $N -anatomy $anatomy -vtk $vtk -vp $vp -dumpfreq $dumpfreq -verbose $verbose
