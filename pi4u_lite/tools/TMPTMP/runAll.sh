#===============================
#  Script for GLIOMA TMCMC
#  * run brain simulaiton
#  * after run likelohood
#===============================
# modules
source /etc/profile.d/modules.sh
module load gcc/4.9
module list

# libs
export MY_BASE=$HOME/GliomaAdvance/lib
export LD_LIBRARY_PATH=$MY_BASE/myVTK/lib/vtk-5.4/:$MY_BASE/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH

# enviroment
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "hi from subtask, we are using $SLURM_CPUS_PER_TASK cores"

# brain simulation set up
program=brain
anatomy=rat
verbose=0
adaptive=1
pID=0
Dscale=20
Tend=11
vtk=0
bDumpIC=0
dumpfreq=1
dumpstart=0
refinefreq=1

# Tumor growth solver
./$program -nthreads $SLURM_CPUS_ON_NODE -anatomy $anatomy -vtk $vtk -dumpstart $dumpstart -dumpfreq $dumpfreq -bDumpIC $bDumpIC -refinefreq $refinefreq -adaptive $adaptive -verbose $verbose -pID $pID -Dscale $Dscale -Tend $Tend


#Likelihood computation
./likelihood

#./cleanLocalTmpDir.sh
