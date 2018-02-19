# set up module, paths
#	TMCMC stuff
export PATH=$PATH:/cluster/home/mavt/chatzidp/usr/mpich3/bin/
export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH

module load gcc

#480 cores -> 10 nodes -> 10x8=80 MPI, 6 threads per MPI rank
#480 cores -> 10 nodes -> 10x12=120 MPI, 4 threads per MPI rank
#336 Cores -> 7 nodes ->  7x12=84 MPI, 4 threads per MPI ranks
cores=480
bsub -J P01V2_4K_J1 -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J P01V2_4K_J2 -w "ended(P01V2_4K_J1)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J P01V2_4K_J3 -w "ended(P01V2_4K_J2)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
