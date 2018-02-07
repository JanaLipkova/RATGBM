#===============================
#  submitting hybrid MPI/OpenMP job
#===============================
#N processors
#M MPI tasks (one per node)
#P openMP threads per node
#-——————————————
N=48 
M=1
P=48
hybrid_program=brain 

# simulation set up
ic=0
vtk=1
ch=0
CHsteps=500
dumpfreq=0.1
#dumpfreq=0.001
tend=2
CFL=0.8


export OMP_NUM_THREADS=$P
mpirun -np $M -pernode ./$hybrid_program -ic $ic -vtk $vtk -nthreads $P -ch $ch -CHsteps $CHsteps -tend $tend -dumpfreq $dumpfreq -CFL $CFL
#bsub -n $N -R ‘span[ptile=$P]’  mpirun -np $M -pernode ./brain -ic 0
