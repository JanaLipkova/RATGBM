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
adaptive=0
verbose=1
ic=1
vtk=1

#parameters
tend=1
dumpfreq=0.1
CFL=0.8
scale=10
kGM=1
kWM=10
kCSF=100


export OMP_NUM_THREADS=$P
bsub -n $N -R 'span[ptile=48]' mpirun -np $M -pernode ./$hybrid_program -nthreads $P -ic $ic -vtk $vtk -tend $tend -dumpfreq $dumpfreq -CFL $CFL -scale $scale -kGM $kGM -kWM $kWM -kCSF $kCSF

