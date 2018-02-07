export OMP_NUM_THREADS=16


for run in 1
do

./submit.sh
#for type in 0 1 2
#do 
#bsub -W 01:00 -n 16 -o out$type ./$run -nthreads 16 -anatomy extended -type $type -adaptive 0 -vtk 1 -tend 0. -dumpfreq 0.1
#mv Brain0000.vtk $run$type.vtk
#done

done
