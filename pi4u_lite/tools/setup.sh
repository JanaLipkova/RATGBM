# Load modules
module purge
module load admin/1.0 lrz/default intel/17.0 mkl/2017 
module load gsl/2.3 blast gcc/4.9
echo "Currently loaded modules:"
module list
# Set up paths
export PATH=$HOME/usr/torc/bin:$PATH
export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH
echo "we use this mpicc:"
which mpicc

# needed on LRZ
export LANG=C
export LC_ALL=C

