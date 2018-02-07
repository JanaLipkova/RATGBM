#===============================
#  submitting OpenMP job
#===============================
#N processors
#P openMP threads per node
#-——————————————
N=48
export OMP_NUM_THREADS=$N
program=brain

# simulation set up
vtk=1
dumpfreq=400
adaptive=0
verbose=1
IC=22

export OMP_NUM_THREADS=$N
#bsub -W 00:10 -n $N ./$program -nthreads $N -vtk $vtk -dumpfreq $dumpfreq -adaptive $adaptive -verbose $verbose  
./$program -nthreads $N -IC $IC -vtk $vtk -dumpfreq $dumpfreq -adaptive $adaptive -verbose $verbose
