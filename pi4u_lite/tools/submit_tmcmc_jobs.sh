#! /bin/bash

# 1) first job - no dependencies
sub_info=$(sbatch run_tmcmc.sh)
echo "submission info: $sub_info"

#get job id from the sub_info
N=4
arr=($sub_info)
IDj1=${arr[N-1]}
echo "j1 id: $IDj1"

#second job, depence: terminatation of job1
sub_info=$(sbatch --dependency=afterany:$IDj1 ./run_tmcmc.sh)
echo $sub_info
arr=($sub_info)
IDj2=${arr[N-1]}
echo "j2: id $IDj2"

#third job, depence: terminatation of job2
#sub_info=$(sbatch --dependency=afterany:$IDj2 ./run_tmcmc.sh)
#echo $sub_info
#arr=($sub_info)
#IDj3=${arr[N-1]}
#echo "j3: id $IDj3"

# check for depenencies of submitted jobs:
squeue --cluster=mpp2 -u $USER -o "%.8A %.4C %.10m %.20E"
