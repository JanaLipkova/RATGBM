# BRUTUS (48 cores per node)
#480 cores -> 10 nodes -> 10x8=80 MPI, 6 threads per MPI rank
#144 cores -> 3 nodes ->  3x12=36 MPI, 4 threads per MPI rank
#336 cores -> 7 nodes ->  7x12=84 MPI, 4 threads per MPI rank
#384 cores -> 8 nodes ->  8x12=96 MPI, 4 threads per MPI rank

#LRZ: (28 cores per node)
#168 cores -> 6 nodes -> 6x7=42 MPI, with 4 threads per MPI rank
MPI=42
ppn=4
mpich_run -n $MPI -ppn $ppn -env TORC_WORKERS 1 -launcher ssh -f $LSB_DJOB_HOSTFILE ./engine_tmcmc



