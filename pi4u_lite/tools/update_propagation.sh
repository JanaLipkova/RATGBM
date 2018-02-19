#!/bin/sh

cases=0
for var in {3866..3870}
do
if [ $(find tmpdir.0.0.0.$var/HGG_data.dat) ]
then
 echo "tmpdir.0.0.0.$var/HGG_data.dat found"
 cases=$((cases+1))
else
 echo "tmpdir.0.0.0.$var/HGG_data.dat not found"
#cd tmpdir.0.0.0.$var/
#echo cp ../tmpdir.0.0.0.3802/runHGG.sh .
#echo sbatch runHGG.sh
#echo $pwd
#cd ..
fi
done

echo "found cases $cases"
