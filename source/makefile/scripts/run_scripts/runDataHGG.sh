#===============================
#  submitting OpenMP job
#===============================
#N processors
#P openMP threads per node
#-——————————————
export LD_LIBRARY_PATH=/cluster/home/mavt/lipkovaj/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/:$LD_LIBRARY_PATH

N=48
export OMP_NUM_THREADS=$N
program=brain

# simulation set up
vtk=1
dumpfreq=400
adaptive=0
verbose=1

export OMP_NUM_THREADS=$N
bsub -W 00:30 -n $N ./$program -nthreads $N -vtk $vtk -dumpfreq $dumpfreq -adaptive $adaptive -verbose $verbose  

