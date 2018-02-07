#=============================
#  running hybrid MPI/openMP allocations
#=============================

Manual to run hybrid MPI/openMP application on
-----------------------------
N processors
M MPI tasks (one per node)
P openMP threads per node
-----------------------------
N = M x P

1) configure your enviroemtn to use P openMP threads:
export OMP_NUM_THREADS=P

2) tell your batch system to give you N processors using P processors per node
bsub -n N -R ‘span[ptile=P]’

3) tell Open MPI to start your program on M nodes using only one MPI task per node
… mpirun -np M -pernode ./hybrid_program

All together:
bsub -n N -R ‘span[ptile=P]’  mpirun -np M -pernode ./hybrid_program


