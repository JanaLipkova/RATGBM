# set up module, paths
#	CMA stuff
module purge
module load gcc/4.7
export PATH=$PATH:/cluster/home/mavt/chatzidp/usr/mpich2/bin/
export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich2/lib/:$LD_LIBRARY_PATH

# submit job
bsub -W 08:00 -n 48 -J GliomaCMA -o GliomaCMA < runcma.lsf.sh

