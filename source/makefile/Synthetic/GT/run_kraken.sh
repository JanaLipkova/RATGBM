Nthreads=8
export OMP_NUM_THREADS=$Nthreads

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
Tend=16


./$program -nthreads $Nthreads -anatomy $anatomy -vtk $vtk -dumpstart $dumpstart -dumpfreq $dumpfreq -bDumpIC $bDumpIC -refinefreq $refinefreq -adaptive $adaptive -verbose $verbose -pID $pID -Dscale $Dscale -Tend $Tend

