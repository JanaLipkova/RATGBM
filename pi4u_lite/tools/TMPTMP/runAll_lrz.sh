
# enviroment
#export OMP_NUM_THREADS=1

# brain simulation set up
program=brain
anatomy=rat
verbose=0
adaptive=1
pID=0
Dscale=10
Tend=11

vtk=0
bDumpIC=0
dumpfreq=1
dumpstart=0
refinefreq=1
ICtype=1

# Tumor growth solver
./$program -nthreads 1 -anatomy $anatomy -vtk $vtk -dumpstart $dumpstart -dumpfreq $dumpfreq -bDumpIC $bDumpIC -refinefreq $refinefreq -adaptive $adaptive -verbose $verbose -pID $pID -Dscale $Dscale -Tend $Tend -ICtype $ICtype


#Likelihood computation
./likelihood

#./cleanLocalTmpDir.sh
