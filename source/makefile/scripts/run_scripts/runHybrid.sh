#============================================
#      Run hybrid MPI/OpenMP job
#============================================

export OMP_NUM_THREADS=2

bsub -n 2 -R 'span[ptile=2]' mpirun -np 1 -pernode ./brain -ic 0


