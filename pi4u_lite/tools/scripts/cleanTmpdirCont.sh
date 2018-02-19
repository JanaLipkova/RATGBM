# remove unneded stuff
for i in tmpdir.0.*
do
  cd $i
  rm brain
  rm HGG_InputParameters.txt
  rm HGG_TumorIC.bin
  rm likelihood
  rm LikelihoodInput.txt
  rm Likelihood.txt
  rm output
  rm read_HGG_TumorIC
  rm README
  rm runAll.sh
  rm testDataCheck.m
  rm Sphere.dat
  rm T1_data.dat
  rm T2_data.dat
  rm PET_data.dat
  cd ../ 
done


