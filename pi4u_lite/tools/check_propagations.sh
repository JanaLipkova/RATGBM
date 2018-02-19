#!/bin/sh

cases=0
for var in {0..4032}
do
if [ $(find tmpdir.0.0.0.$var/HGG_data.dat) ]
then
#echo "tmpdir.0.0.0.$var/HGG_data.dat found"
cases=$((cases+1))
#echo $cases
else
echo "tmpdir.0.0.0.$var/HGG_data.dat not found"
fi

done

echo "found cases $cases"
