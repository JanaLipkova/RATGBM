#288 cores -> 6 nodes -> 6x8=48 MPI, 6 threads per MPI rank
#480 cores -> 10 nodes -> 10x8=80 MPI, 6 threads per MPI rank
#144 cores -> 3 nodes  ->  3x12=36 MPI, 4 threads per MPI rank
#288 core  -> 6 nodes  ->  6x12=48 MPI, 4 threads per MPI rank 
#336 cores -> 7 nodes ->  7x12=84 MPI, 4 threads per MPI rank
MPI=84
ppn=12
mpich_run -n $MPI -ppn $ppn -env TORC_WORKERS 1 -launcher ssh -f $LSB_DJOB_HOSTFILE ./propagation_tool curgen_db_022.txt



