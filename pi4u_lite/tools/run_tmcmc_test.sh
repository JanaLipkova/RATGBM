#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J SBall
#SBATCH --clusters=mpp1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=2
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=00:10:00

source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default intel/16.0 mkl/11.3
module load gsl/2.3 blast gcc/4.8
module list

export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH
echo "we use this mpicc:"
which mpicc

export PATH=$HOME/usr/torc/bin:$PATH
export LANG=C
export LC_ALL=C
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo "In the directory: $PWD"
echo "Running program on $SLURM_CPUS_ON_NODE nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."

#LRZ: (28 cores per node)
#  1 cores ->  28 nodes  ->   1x14=14 MPI, each with 2 threads per MPI ran
#  5 cores -> 140 nodes  ->   5x14=70 MPI, each with 2 threasa
# 10 cores -> 280 nodes  -> 10x14=140 MPI, each with 2 threads  
ppn=8
mpirun -np $SLURM_NTASKS -ppn $ppn -env TORC_WORKERS 1 ./engine_tmcmc

#mpich_run -n $MPI -ppn $ppn -env TORC_WORKERS 1 -launcher ssh -f $LSB_DJOB_HOSTFILE ./engine_tmcmc



#bsub -J SNAIC_4K_J3 -w "ended(SNAIC_4K_J2)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh
