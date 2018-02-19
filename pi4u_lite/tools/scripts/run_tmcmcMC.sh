# set up module, paths
#	TMCMC stuff
export PATH=$PATH:/cluster/home/mavt/chatzidp/usr/mpich3/bin/
export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
export PATH=$HOME/usr/torc/bin:$PATH

module load gcc

#480 cores -> 10 nodes -> 10x8=80 MPI, 6 threads per MPI rank
#144 cores -> 3 nodes ->  3x12=36 MPI, 4 threads per MPI rank 
cores=336
bsub -J SNAIC_4K_J1 -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J SNAIC_4K_J2 -w "ended(SNAIC_4K_J1)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
bsub -J SNAIC_4K_J3 -w "ended(SNAIC_4K_J2)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
#bsub -J SNAIC_4K_J4 -w "ended(SNAIC_4K_J3)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
