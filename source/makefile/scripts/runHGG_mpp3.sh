#!/bin/bash
#SBATCH -o myjob.%j.%N.out
#SBATCH -D .
#SBATCH -J P00
#SBATCH --clusters=mpp3
#SBATCH --nodes=1-1
#SBATCH --cpus-per-task=48
# if mpp1 is used replace above 28 by 16
#SBATCH --export=NONE
#SBATCH --time=00:20:00

source /etc/profile.d/modules.sh
module load gcc/4.9

export MY_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$MY_BASE/myVTK/lib/vtk-5.4/:$MY_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# brain simulation set up
program=brain
anatomy=rat
verbose=1
adaptive=1
pID=0
Dscale=20

vtk=1
bDumpIC=0
dumpfreq=1
dumpstart=0
refinefreq=1
UQtype=0

echo "In the directory: $PWD"
echo "Running program on $SLURM_CPUS_ON_NODE nodes, each with $SLURM_CPUS_PER_TASK cores."

./$program -nthreads $SLURM_CPUS_ON_NODE -anatomy $anatomy -vtk $vtk -dumpstart $dumpstart -dumpfreq $dumpfreq -bDumpIC $bDumpIC -refinefreq $refinefreq -adaptive $adaptive -verbose $verbose -pID $pID -UQtype $UQtype -Dscale $Dscale 

