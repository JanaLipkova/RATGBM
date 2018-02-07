#==================================
#    Cahn Hiliard solver 
#   compute phase filed function
#==================================
P=48
vtk=1
steps=400
width=3
verbose=1
dumpfreq=50

export OMP_NUM_THREADS=$P
bsub -W 02:00 -n $P ./brain -nthreads $P -vtk $vtk -steps $steps -width $width -verbose $verbose -dumpfreq $dumpfreq 
