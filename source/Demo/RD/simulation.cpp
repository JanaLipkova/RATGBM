/*
 *  simulation.cpp
 *  RD
 *
 *  Created by Diego Rossinelli on 1/13/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "MRAGHeaders.h"

using namespace MRAG;

#include "simulation.h"

struct CCell
{
	float u[2];
	float dudt[2];
	float tmp[2];
	
	CCell()
	{
		u[0] = u[1] = 0;
		dudt[0] = dudt[1] = 0;
		tmp[0] = tmp[1] = 0;
	}
	
	CCell(float){abort();}
	
	CCell(float u_, float v_)
	{
		u[0] = u_;
		u[1] = v_;
		dudt[0] = dudt[1] = 0;
		tmp[0] = tmp[1] = 0;
	}
	
	void operator += (CCell t) 
	{
		for(int i=0; i<2; i++) u[i] += t.u[i];
		for(int i=0; i<2; i++) dudt[i] += t.dudt[i];
	}
	
	operator Real() 
	{
		return (Real)u[0];
	}
	
	void integrate(float dt)
	{
		for(int i=0; i<2; i++) u[i] += dt*dudt[i];
		for(int i=0; i<2; i++) dudt[i] = tmp[i] = 0;
	}
	
	template<int i>
	Real evaluate_concentration(double dt)
	{
		return  u[i]+dt*dudt[i];
	}
};

CCell operator*(const CCell& p, Real v)
{
	CCell c;
	for(int i=0; i<2; i++) c.u[i] = p.u[i]*v;
	for(int i=0; i<2; i++) c.dudt[i] = p.dudt[i]*v;
	return c;
}

#ifdef _MRAG_GLUT_VIZ
RGBA convertToRGBA(CCell& p)
{
	const double c = p.u[1]; 
	const double R = max(0., min(1., 2*(c-0.5)));
	const double G = max(0., min(1., 1-2*fabs(c-0.5)));
	const double B = max(0., min(1., 1-2*c));
	RGBA color(R,G,B,0);
	return color;
}
#endif

template <typename T, int i> inline Real RD_projector_impl(const T&t)
{
	return (Real)(t.u[i]);
}

make_projector(RD_projector, RD_projector_impl)

static const int resJump = 1;
const double refinement_tolerance = 1e-3;
const double compression_tolerance = 1e-4;
const double grayscott_F = 0.04;
const double grayscott_kappa = 0.06;
const double grayscott_diffusion_rates[2] = {2.e-06, 1.e-06};

typedef Block< CCell, blockSize, blockSize, blockSize> B;
typedef _WAVELET_TYPE W;
#ifdef _MRAG_TBB
    typedef Multithreading::BlockProcessing_Pipeline_TBB<B, BlockLab, nThreads+1> BlockProcessing;
#else
    typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
#endif

struct Integrate_TR
{
	double dt;
	
	Integrate_TR(const double dt_): dt(dt_){}
	Integrate_TR(const Integrate_TR& c): dt(c.dt) {}
	
	template<typename B>
	inline void operator()(const BlockInfo& info, B& b) const
	{
		typedef typename B::ElementType E;
		
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		
		E* ptrE = &(b[0]);
		
		for(int iE=0; iE<n; iE++, ptrE++)
			ptrE->integrate(dt);
	}
};

struct UpdateRHS_TR
{
	double dt;
	
	UpdateRHS_TR(const double dt_): dt(dt_){}
	UpdateRHS_TR(const UpdateRHS_TR& c): dt(c.dt) {}
	
	template<typename B>
	inline void operator()(const BlockInfo& info, B& b) const
	{
		typedef typename B::ElementType E;
		
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		
		E* ptrE = &(b[0]);
		
		for(int iE=0; iE<n; iE++, ptrE++)
		{
			for(int i=0;i<2;i++) ptrE->dudt[i] = ptrE->tmp[i];
			for(int i=0;i<2;i++) ptrE->u[i] -= dt*ptrE->tmp[i];
		}		
	}
};

struct ComputeRHS_TR
{
	int stencil_start[3], stencil_end[3];
	Real dt;
	
	ComputeRHS_TR(const ComputeRHS_TR& c): dt(0.)
	{
		memcpy(this, &c, sizeof(ComputeRHS_TR) );
	}
	
	ComputeRHS_TR(Real dt_): dt(dt_)
	{
		stencil_start[0] = stencil_start[1] = -1; stencil_start[2] = 0;
		stencil_end[0] = stencil_end[1] = 2; stencil_end[2] = 1;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
	{	
		typedef BlockType B;
		typedef typename BlockType::ElementType E;
		
		const Real h = info.h[0];
		const Real lapl_factor = 1/(h*h);
		const Real nu_U = grayscott_diffusion_rates[0];
		const Real nu_V = grayscott_diffusion_rates[1];
		const Real ampl_r = 2.0;
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
				{	
					const Real currU = i(ix, iy, iz).template evaluate_concentration<0>(dt);
					const Real currV = i(ix, iy, iz).template evaluate_concentration<1>(dt);
					
					o(ix,iy,iz).tmp[0] = nu_U*lapl_factor*(
						i(ix+1, iy, iz).template evaluate_concentration<0>(dt)+
						i(ix-1, iy, iz).template evaluate_concentration<0>(dt)+
						i(ix, iy+1, iz).template evaluate_concentration<0>(dt)+
						i(ix, iy-1, iz).template evaluate_concentration<0>(dt)+
						i(ix, iy, iz+1).template evaluate_concentration<0>(dt)+
						i(ix, iy, iz-1).template evaluate_concentration<0>(dt) - 6*currU)
						+
						ampl_r*(-currU*currV*currV + grayscott_F*(1 - currU));
					
					o(ix,iy,iz).tmp[1] = nu_V*lapl_factor*(
						i(ix+1, iy, iz).template evaluate_concentration<1>(dt)+
						i(ix-1, iy, iz).template evaluate_concentration<1>(dt)+
						i(ix, iy+1, iz).template evaluate_concentration<1>(dt)+
						i(ix, iy-1, iz).template evaluate_concentration<1>(dt)+
						i(ix, iy, iz+1).template evaluate_concentration<1>(dt)+
						i(ix, iy, iz-1).template evaluate_concentration<1>(dt) - 6*currV)
						+
						ampl_r*(currU*currV*currV - (grayscott_F + grayscott_kappa)*currV);
				}
	}
};

Grid<W, B> * grid;
BlockProcessing blockProcessing;
Refiner_SpaceExtension refiner(resJump);
Compressor compressor(resJump);
BlockFWT<W, B, RD_projector> blockfwt;
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
void _ic(Grid<W, B>& grid)
{
	vector<BlockInfo> vInfo = grid.getBlocksInfo();

	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo& info = vInfo[i];
		B& block = grid.getBlockCollection()[info.blockID];
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
				{
					double x[3];
					info.pos(x, ix, iy,iz);
					
					const double r[3] = {
						x[0] - 0.5, x[1]-0.5, x[2]-0.5
					};
					const double IrI = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
					const double center =0.05;
					const double hs = 0.01;
					const double eps = hs*2;
					double fInside = _smartW(center-hs, eps, IrI, 1)*_smartW(center+hs, eps, IrI, -1);
					
					fInside = pow(fInside,4);
					const double noise_factor = 0.1;
					block(ix,iy,iz).u[0] = fInside*(0.5+0.005*noise_factor*drand48())  + (1-fInside)*(1.0 + 0.01*noise_factor*drand48()); 
					block(ix,iy,iz).u[1] = fInside*(0.25+0.0025*noise_factor*drand48()) + (1-fInside)*(0.0 + 0.01*noise_factor*drand48()); 
					block(ix,iy,iz).dudt[0] = 0;
					block(ix,iy,iz).dudt[1] = 0;
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
	
	printf("RD INIT! nThreads=%d, blockSize=%d Wavelets=w%s (blocksPerDimension=%d, maxLevel=%d)\n",
		   nThreads, blockSize, wavelet_name.data(), blocksPerDimension, maxLevel);
	Environment::setup();
	
	grid = new Grid<W,B>(blocksPerDimension,blocksPerDimension, blocksPerDimension, maxStencil);
	grid->setCompressor(&compressor);
	grid->setRefiner(&refiner);
	stSorter.connect(*grid);
	
	_ic(*grid);
	Science::AutomaticRefinement<0,0>(*grid, blockfwt, refinement_tolerance, maxLevel, 1, &profiler, _ic);
	_ic(*grid);
	Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, 1, &profiler, _ic);
}

void _step(float t, float dt, BoundaryInfo* boundaryInfo, const int nParallelGranularity)
{
	//stSorter.startSession(dt, 4, 0);
	ComputeRHS_TR task1(dt);
	UpdateRHS_TR task2(dt);
	Integrate_TR task3(dt);

	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const BlockCollection<B>& collecton =  grid->getBlockCollection();
	
	profiler.getAgent("computing").start();
	for(int i=0; i<30; i++)
	{
		blockProcessing.pipeline_process(vInfo, collecton, *boundaryInfo, task1);
		BlockProcessing::process(vInfo, collecton, task2, nParallelGranularity);
		BlockProcessing::process(vInfo, collecton, task3, nParallelGranularity);
	}
	profiler.getAgent("computing").stop();
	/*while(true)
	{
		double currTime, currDeltaT;
		int level;
		SpaceTimeSorter::ETimeInterval type;
		vector<BlockInfo> vInfo;
		
		const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
		printf("block to process in parallel %d\n", vInfo.size());
		if (type == SpaceTimeSorter::ETimeInterval_Start)
		{
			ComputeRHS_TR task1(currTime);
			UpdateRHS_TR task2(currTime);
			
			profiler.getAgent("computing").start();
			blockProcessing.pipeline_process(vInfo, grid->getBlockCollection(), *boundaryInfo, task1);
			BlockProcessing::process(vInfo, grid->getBlockCollection(), task2, nParallelGranularity);
			profiler.getAgent("computing").stop();
		}
		else
		{
			Integrate_TR task(currTime);
			profiler.getAgent("computing").start();
			BlockProcessing::process(vInfo, grid->getBlockCollection(), task, nParallelGranularity);
			profiler.getAgent("computing").stop();
		}
		
		if (!bContinue) break;
	}
	stSorter.endSession();		*/
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
	dt = 5.0*dx*dx*(32*32);
	
	printf("STEP:%02d computing phase now\n", static_counter);
	_step(t,dt, boundaryInfo, nParallelGranularity);
	_step(t,dt, boundaryInfo, nParallelGranularity);
	t+=2*dt;
	
	profiler.printSummary();
	
	//Science::AutomaticCompression<0,0>(*grid, blockfwt, compression_tolerance, -1, &profiler);
	printMemorySummary(static_counter, "compressing after computing");
	static_counter++;
}

#ifdef _MRAG_GLUT_VIZ
GridViewer viewer(!W::bIsCellCentered, true);
#endif
void simulation_render(bool bDrawTextures)
{
#ifdef _MRAG_GLUT_VIZ
	viewer.drawContent(*grid, grid->getBlockCollection());
	viewer.drawSketch(*grid,false);
#endif
}
