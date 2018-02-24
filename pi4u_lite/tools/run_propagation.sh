#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J P16p
#SBATCH --clusters=mpp2
#SBATCH --ntasks=224
#SBATCH --cpus-per-task=2
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=25:25:00

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
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export MY_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$MY_BASE/myVTK/lib/vtk-5.4/:$MY_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH

echo "In the directory: $PWD"
echo "Running program with $SLURM_NTASKS MPI tasks, each with $SLURM_CPUS_PER_TASK cores."

#echo "Running program on $SLURM_CPUS_ON_NODE nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."

#----------------------------
#LRZ: MPP2 (28 cores per node)
#---------------------------
#  1 cores ->  28 nodes  ->   1x14=14 MPI, each with 2 threads per MPI ran
# 10 cores -> 280 nodes  -> 10x14=140 MPI, each with 2 threads  
# ppn=14
# -------------------------------
#LRZ: MPP1 (16 cores per node)
# if LR 128^3 grid: -> 8 MPI per node:
# if HR: 256^3 grid -> 4 MPI per node: - otherwise doesn't fit in the memory
#-----------------------------------
# with 2 threads per MPI -> 16/2 = 8 MPI per node, i.e. ppn=8
#  1 cores ->  16 nodes  ->   1x8=8 MPI, each with 2 threads per MPI ran
# 10 cores -> 160 nodes  -> 10x8=80 MPI, each with 2 threads
# ppn=8
#-----------------------------
# with 4 threads per MPI -> 16/4 = 4 MPI per node, i.e. ppn=4
#  1 cores ->  16 nodes  ->   1x4=4 MPI, each with 4 threads per MPI ran
# 10 cores -> 160 nodes  -> 10x8=80 MPI, each with 4 threads
#-----------------------------
ppn=14
mpirun -np $SLURM_NTASKS -ppn $ppn -env TORC_WORKERS 1 ./propagation_tool curgen_db_022.txt 


#bsub -J SNAIC_4K_J3 -w "ended(SNAIC_4K_J2)" -W 10:00 -n $cores -R "span[ptile=48]" -R "select[model==Opteron6174]" < runtmcmc.lsf.sh