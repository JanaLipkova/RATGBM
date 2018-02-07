#===============================
#  Script for GLIOMA TMCMC
#  * run brain simulaiton
#  * after run likelohood
#===============================

N=8
export OMP_NUM_THREADS=$N

# brain simulation set up
program=brain
vtk=0
dumpfreq=50
adaptive=1
verbose=0

./$program -nthreads $N -vtk $vtk -dumpfreq $dumpfreq -adaptive $adaptive -verbose $verbose  

#likelihood
./likelihood

