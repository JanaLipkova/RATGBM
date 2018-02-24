echo "we use this mpicc:"
which mpicc

export PATH=$HOME/usr/torc/bin:$PATH

mpirun -np 32 ./engine_tmcmc

#mpich_run -n $MPI -ppn $ppn -env TORC_WORKERS 1 -launcher ssh -f $LSB_DJOB_HOSTFILE ./engine_tmcmc

