#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J M00_5
#SBATCH --clusters=mpp3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=06:00:00

# Modules
source /etc/profile.d/modules.sh
module purge
module load admin/1.0 lrz/default 
module load gsl/2.3 gcc/4.9
module list

# Libraries
export PATH=$PATH:$HOME/pi4u-libs/mpich-install/bin/
export LD_LIBRARY_PATH=$HOME/pi4u-libs/mpich-install/lib/:$LD_LIBRARY_PATH
echo "we use this mpicc:"
which mpicc

export PATH=$HOME/usr/torc/bin:$PATH
export LD_LIBRARY_PATH=$HOME/usr/torc/bin:$LD_LIBRARY_PATH

export MY_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$MY_BASE/myVTK/lib/vtk-5.4/:$MY_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/hpc/txh01/di49zin/ssm-libs/boost_1_63_0/lib/:$LD_LIBRARY_PATH


export LANG=C
export LC_ALL=C
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE  

echo "In the directory: $PWD"
echo "Running program on $SLURM_TASKS nodes with $SLURM_CPUS_PER_TASK tasks, each with $SLURM_CPUS_PER_TASK cores."

#----------------------------
#LRZ: MPP2 (28 cores per node)
#---------------------------
#  1 cores ->  28 nodes  ->   1x14=14 MPI, each with 2 threads per MPI ran
# 10 cores -> 280 nodes  -> 10x14=140 MPI, each with 2 threads
# ppn=14
# -------------------------------

mpirun -np 64 -env TORC_WORKERS 1 ./engine_tmcmc

#mpich_run -n $MPI -ppn $ppn -env TORC_WORKERS 1 -launcher ssh -f $LSB_DJOB_HOSTFILE ./engine_tmcmc

