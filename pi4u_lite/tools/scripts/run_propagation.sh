# set up module, paths
#	TMCMC stuff
export PATH=$PATH:/cluster/home/mavt/chatzidp/usr/mpich3/bin/
export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH

module load gcc

#288 cores -> 6 nodes -> 6x8=48 MPI, 6 threads per MPI rank
#480 cores -> 10 nodes -> 10x8=80 MPI, 6 threads per MPI rank
#144 cores -> 3 nodes ->  3x12=36 MPI, 4 threads per MPI rank
#288 core  -> 6 nodes  ->  6x12=48 MPI, 4 threads per MPI rank
#336 cores -> 7 nodes ->  7x12=84 MPI, 4 threads per MPI rank
cores=336
bsub -J P01_prop4K -w "ended(P01_propHR_4K)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runpropagation.lsf.sh


