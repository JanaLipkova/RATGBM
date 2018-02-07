==========================================
	    GLIOMA SIMULATION 
		SDE3.h
==========================================


==========================================
  	BRUTUS (Linux):
==========================================

MODULES:

Glioma code work with gcc/4.4.6 not the 4.7.2 thus:

module unload gcc/4.7.2
module load gcc/4.4.6

export LD_LIBRARY_PATH=/cluster/home/mavt/lipkovaj/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/work/infk/cconti/hypre-2.10.0b/src/hypre/lib/:$LD_LIBRARY_PATH

==========================================

==========================================
   on Jana's MAC:
==========================================
export DYLD_LIBRARY_PATH=/opt/intel/Compiler/11.1/084/Frameworks/tbb/ia32/cc4.0.1_os10.5.4/lib/:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/usr/local/lib/vtk-5.6/:$DYLD_LIBRARY_PATH
==========================================


COMPILATAION

 > make clean
 > make -j 8

   or to avoid assertion checks:
> make config=production -j 144
	
 >> creates executable brain

-----------------------------------------
	Glioma_extendedAnatomy.cpp
-----------------------------------------
run:

./brain -nthreads 2 -anatomy extended -type 2 -adaptive 0 -vtk 1 -tend 0.1 -dumpfreq 0.1

./brain -nthreads 2 -anatomy advection -adaptive 0 -vtk 1 -tend 0.1 -dumpfreq 10 -dt 0.01 -CFL 0.5 -ic 1 -RK2 1

./brain -nthreads 2 -anatomy laplace -adaptive 0 -tol 1e-3 -dumpfreq 5 -vtk 1 -ic 0

./brain -nthreads 2 -anatomy pressure -adaptive 0 -tol 1e-3 -dumpfreq 5 -vtk 1 


-D Subject04       = define to use anatomy of subject04 from BrainWeb (skull + bone marrow, wm + gm, glioma) or OLD_VERSION
-nthreads          = #threads for parallelisation
-anatomy           = to choce between basic (wm+gm), or extended (skull+bone marrow)
-vtk 	           = 1 to enable dumping data to vtk format, 0 without vtk output
-type              = 0 original BrainWeb, 1 Subject04 (fuzzy), 2 Subject04 (discrete)
-adaptive          = 1 for adaptive grid, 0 for uniform grid
-tend              = final time, [days]
-dumpfreq          = how often to dump data, also how often to refine


-----------------------------------------
	Original version RD/SDE3D
-----------------------------------------
EXECUTION PARAMERERS
./brain -nthreads 48 -ic 0 -vtk 1 -tend 54 -dumpfreq 0.1 

-nthreads 	= #threads for parallelisation
-ic       	= type of initial condition, ic = 0 for BrainWeb anatomy, ic = 1 for patient specific anatomy
-vtk	 	= 1 to enable dumping data to vtk format, 0 without vtk output
-tend  	 	= final time, [days]
-dumpfreq 	= how often to dump data, also how often to refine
-unifrom 	= 1 to use uniform grid (needed e.g. for UQ), 0 to use adptive grid(refinemnet and compression)
-binary		= 1 to dump binary output, 0 no	

==========================================
    Notes to compilation
==========================================

= on MAC one needs DYLD_LIBRARY_PATH insread of LD_LIBRARY_PATH 
= to check that path was expoerted correctly do:
	echo $DYLD_LIBRARY_PATH

==========================================


