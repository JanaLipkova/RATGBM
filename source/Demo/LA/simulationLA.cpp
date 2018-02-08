/*
 *  simulation.cpp
 *  RD
 *
 *  Created by Diego Rossinelli on 1/13/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "LA_MRAGHeaders.h"

using namespace MRAG;

#include "simulationLA.h"

struct PLS
{
	float rho,rho_0;
	float drho_dt;
	float tmp;

	PLS(): rho(0), drho_dt(0), tmp(0), rho_0(0) {}
	PLS(float rho_, float drho_dt_): rho(rho_), drho_dt(drho_dt_),tmp(0), rho_0(0) {}
	PLS(float rho_): rho(rho_), drho_dt(0),tmp(0), rho_0(0) {}
	
	Real levelset() const{return rho;}
	
	void operator += (PLS t) 
	{
		rho_0 +=t.rho_0;
		rho += t.rho;
		drho_dt += t.drho_dt;
	}
	
	operator Real() {return (Real)rho;}
	
	void integrate(float dt)
	{
		rho += drho_dt*dt;
		drho_dt = 0;
		tmp = 0;
	}
	
	float evaluate_rho(double dt){return rho + drho_dt*dt;}
};

PLS operator*(const PLS& p, Real v)
{
	PLS t(p.rho*v, p.drho_dt*v);
	return t;
}

template <typename T, int i> inline Real levelset_projector_impl(const T&t)
{
	return (Real)t.rho;
}

make_projector(LA_projector, SimpleLevelsetBlock<PLS>::levelset_projector_impl)


template<typename RealType>
static void _velocity(const RealType x[3],  RealType  t, RealType  v[3])
{
	const double factor = 1.0;
	
	v[0] =  factor*2.0* pow(sin(M_PI*x[0]), 2) * sin(2*M_PI*x[1]) * sin(2*M_PI*x[2]);
	v[1] = -factor* pow(sin(M_PI*x[1]), 2) * sin(2*M_PI*x[0]) *  sin(2*M_PI*x[2]);
	v[2] = -factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1]);

}

#ifdef _MRAG_GLUT_VIZ
RGBA convertToRGBA(PLS& p)
{
	const double c = p.rho; 
	const double R = max(0., min(1., 2*(c-0.5)));
	const double G = max(0., min(1., 1-2*fabs(c-0.5)));
	const double B = max(0., min(1., 1-2*c));
	RGBA color(R,G,B,0);
	return color;
}
#endif

#include "LA_blockprocessings.h"

static const int resJump = 1;
static const int narrowBandWidth = 6;
		
typedef SimpleLevelsetBlock< PLS, narrowBandWidth, blockSize, blockSize, blockSize> B;
typedef _WAVELET_TYPE W;
typedef Multithreading::BlockProcessing_Pipeline_TBB<B, BlockLab, nThreads+2> BlockProcessing;

Grid<W, B> * grid;
BlockProcessing blockProcessing;
Refiner_SpaceExtension refiner(resJump);
Compressor compressor(resJump);
BlockFWT<W, B, LA_projector> blockfwt;
SpaceTimeSorter stSorter;
Profiler profiler;

inline const double _smartW(double x0, double eps, double x, int dir=1)
{
	const double alpha = M_PI*min(1., max(0., (x-x0+eps*0.5)/eps));
	assert(alpha>=0 && alpha<=M_PI);
	if (dir == 1) return 0.5+0.5*sin(alpha-M_PI/2);
	else return 0.5+0.5*cos(alpha);
}
#ifndef drand48
#define drand48() (rand()/(double)(RAND_MAX))
#endif

static float _ic_func4(float x[3]) 
{
	const float r1 = sqrt(pow(x[0]-0.35,2) + pow(x[1]-0.35, 2) + pow(x[2]-0.35,2));
	const double d1= 0.15-r1;

	const float r2 = sqrt(pow(1-x[0]-0.35,2) + pow(1-x[1]-0.35, 2) + pow(1-x[2]-0.35,2));
	const double d2= 0.15-r2;

	return fabs(d1)<fabs(d2)?d1:d2;
}

static void _ic(Grid<W, B>& grid) 
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo& info = vInfo[i];
		B& block = grid.getBlockCollection().lock(info.blockID);
		
		block.setH(info.h[0]);
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
				{
					float x[3];
					
					info.pos(x, ix, iy, iz);
					
					block(ix,iy,iz).rho = _ic_func4(x);
					block(ix,iy,iz).rho_0 = 0;
					block(ix,iy,iz).tmp = 0;
				}
		
		grid.getBlockCollection().release(info.blockID);
	}
}

void simulation_init()
{
	std::map<string, string> mapNames;
	mapNames[string(typeid(Wavelets_AverageInterp5thOrder).name())] = string("5");
	mapNames[string(typeid(Wavelets_Interp4thOrder).name())] = string("4");
	mapNames[string(typeid(Wavelets_AverageInterp3rdOrder).name())] = string("3");
	mapNames[string(typeid(Wavelets_Interp2ndOrder).name())] = string("2");
	const string wavelet_name = mapNames[string(typeid(W).name())];
	
	printf("suggested commands:\n");
	printf("mv test test_t%d_b%d_w%s\n", nThreads, blockSize, wavelet_name.data());
	printf("bsub -n %d -W 00:40 -o output_test_t%d_b%d_w%s ./test_t%d_b%d_w%s\n", 
		   nThreads, 
		   nThreads, blockSize, wavelet_name.data(),
		   nThreads, blockSize, wavelet_name.data());
	
	printf("LA INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n",
		   nThreads, blockSize, wavelet_name.data(), blocksPerDimension, maxLevel);
	Environment::setup();
	
	grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
	grid->setCompressor(&compressor);
	grid->setRefiner(&refiner);
	stSorter.connect(*grid);

	profiler.getAgent("initial condition").start();
	_ic(*grid);
	Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _ic);
	_ic(*grid);
	Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, 1, &profiler, _ic);
	profiler.getAgent("initial condition").stop();
}

void _step(float t, float dt, BoundaryInfo* boundaryInfo, const int nParallelGranularity)
{
	static const bool bUsePipeline = false;
	
	vector<BlockInfo> vInfo =grid->getBlocksInfo(); 

	{
		LA_ComputeRHS_TR_WENO5 task1(dt, 1.);
		LA_UpdateRHS_TR task2(0);
		LA_Integrate_TR task3(dt);	
		
		profiler.getAgent("computing").start();
		{	
			if (bUsePipeline)
				blockProcessing.pipeline_process(vInfo, grid->getBlockCollection(), *boundaryInfo, task1);
			else
				blockProcessing.process(vInfo, grid->getBlockCollection(), *boundaryInfo, task1,1);
			
			blockProcessing.process(vInfo, grid->getBlockCollection(), task2,nParallelGranularity);
			blockProcessing.process(vInfo, grid->getBlockCollection(), task3,nParallelGranularity);
		}
		profiler.getAgent("computing").stop();
	}

	{
		LA_ComputeRHS_WENO5_REINIT computeRhs;
		LA_Integrate_REINIT integrate(1e-4);
		
		profiler.getAgent("computing").start();
		for(int i=0; i<5; i++)
		{
			if (bUsePipeline)
				blockProcessing.pipeline_process(vInfo, grid->getBlockCollection(), *boundaryInfo, computeRhs);
			else
				blockProcessing.process(vInfo, grid->getBlockCollection(), *boundaryInfo, computeRhs,1);
			
			blockProcessing.process(vInfo, grid->getBlockCollection(),  integrate,nParallelGranularity);
		}
		profiler.getAgent("computing").stop();
	}
}

void printMemorySummary(int step_id, const char * stage)
{
	printf("Memory report for step (%02d), stage is=%s\n", step_id, stage);

	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	printf("Total Memory: %.2f MB, Boundaries: %.2f MB #blocks = %d\n", 
		   grid->getMemorySize(), grid->getBoundaryInfo().getMemorySize(), vInfo.size());
}

void simulation_run(int nsteps)
{
	static int static_counter = 0;
	
	//Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler);
	printMemorySummary(static_counter, "refining before computing");
	
	const int nParallelGranularity = (grid->getBlocksInfo().size()<=8 ? 1 : 4);
	BoundaryInfo* boundaryInfo = &grid->getBoundaryInfo();
	double dx = 1./(blockSize*pow(2., maxLevel));
	const double CFL = 0.45;
	dt = 2.0*dx*CFL;
	
	printf("STEP:%02d computing phase now\n", static_counter);
	_step(t,dt, boundaryInfo, nParallelGranularity);
	t+=dt;
	
	profiler.printSummary();
	
	//Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
	printMemorySummary(static_counter, "compressing after computing");
	static_counter++;
}

#ifdef _MRAG_GLUT_VIZ
GridViewer viewer(!W::bIsCellCentered, true);
void _drawLevelsetIntersections()
{
	const double iso_val = 0.0;
	glColor3f(1.0,0.0,0.573);

	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo& info = vInfo[i];
		const Real h = info.h[0];

		B& block = grid->getBlockCollection().lock(info.blockID);
		
		glPointSize(300*h);

		glBegin(GL_POINTS);
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY-1; iy++)
				for(int ix=0; ix<B::sizeX-1; ix++)
				{
					const double uC = block(ix,iy,iz).rho-iso_val;
					const double uR = block(ix+1,iy,iz).rho-iso_val;
					const double uT = block(ix,iy+1,iz).rho-iso_val;

					const bool vCrossFlag[2] = {uC*uR<0., uC*uT<0.};

					if (vCrossFlag[0] || vCrossFlag[1])
					{
						const double omega[2] = {
							vCrossFlag[0]? uC/fabs(uC-uR) : 0,
							vCrossFlag[1]? uC/fabs(uC-uT) : 0
						};

						double x[3] = {0,0,0};
						info.pos(x, ix, iy, iz);

						glVertex2f(x[0]+omega[0]*h, x[1]+omega[1]*h);
					}
				}
		glEnd();

		grid->getBlockCollection().release(info.blockID);
	}
}
#endif
void simulation_render(bool bDrawTextures)
{
#ifdef _MRAG_GLUT_VIZ
	_drawLevelsetIntersections();
	//viewer.drawContent(*grid, grid->getBlockCollection());
	viewer.drawSketch(*grid,false);
#endif
}
