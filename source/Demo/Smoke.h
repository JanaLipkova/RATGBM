/*
 *  Smoke.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/5/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <string>


#include "MRAGcore/MRAGCommon.h"
#include "MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp3rdOrder.h"
#include "MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "MRAGcore/MRAGWavelets_AverageInterp5thOrder.h"
#include "MRAGcore/MRAGBlock.h"
#include "MRAGcore/MRAGrid.h"
#include "MRAGcore/MRAGBlockFWT.h"
#include "../MRAGscience/MRAGScienceCore.h"
#include "../MRAGscience/MRAGSpaceTimeSorter.h"
#include "../MRAGscience/MRAGRefiner_SpaceExtension.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_TBB.h"
#include "MRAGcore/MRAGProfiler.h"
#include "../MRAGvisual/MRAGVisualTypes.h"



using namespace MRAG;

struct P
{
	float rho;
	float drho_dt;
	float tmp;
	
	P(): rho(0), drho_dt(0) , tmp(0) {}
	
	P(float rho_, float drho_dt_): rho(rho_), drho_dt(drho_dt_), tmp() {}
	
	void operator += (P t)
	{
		rho += t.rho;
		drho_dt += t.drho_dt;
	}
	
	operator Real()
	{
		return (Real)rho;
	}
	
	void integrate(float dt)
	{
		rho += drho_dt*dt;
		drho_dt = 0;
	}
	
	float evaluate_rho(double dt)
	{
		return rho + drho_dt*dt;
	}
};



RGBA convertToRGBA(P& p)
{
	RGBA c(max(0.0, 2.0*p.rho - 0.5), 1.0 - 2.0*fabs(p.rho - 0.5),max(0.0, 1.0-2.0*p.rho),0) ;
	
	c.r *= p.rho*2;
	c.g *= p.rho*2;
	c.b *= p.rho*2;
	
	return c;
}


P operator*(const P& p, Real v)
{
	P t(p.rho*v, p.drho_dt*v);
	return t;
}

template <typename T, int i> inline Real smoke_projector_impl(const T&t)
{
	return (Real)t.rho;
}

make_projector(smoke_projector, smoke_projector_impl)


template <int nDim, typename T>
void _velocity(T x[3], T t, T v[3]) 
{
	
	if (nDim==3)
	{
		const double factor = 0.2;
		v[0] =  factor*2.0* pow(sin(M_PI*x[0]), 2) * sin(2*M_PI*x[1]) * sin(2*M_PI*x[2]);
		v[1] = -factor* pow(sin(M_PI*x[1]), 2) * sin(2*M_PI*x[0]) *  sin(2*M_PI*x[2]);
		v[2] = -factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1]);
	}
	else
	{
		const double factor = t>10.0? -0.2:0.2;
		const double p[3] = { 2*M_PI*(x[0] - 0.00), 2*M_PI*(x[1] - 0.00), 2*M_PI*(x[2] - 0.00)};
		v[0] = -factor*2.0* pow(sin(p[0]*0.5), 2) * sin(p[1]*0.5) *cos(p[1]*0.5);
		v[1] = 2.0*factor* pow(sin(p[1]*0.5), 2) * sin(p[0]*0.5) *cos(p[0]*0.5);
		v[2] =  0;
	}
}


template<int nDim=3>
struct ComputeRHS
{
	int stencil_start[3], stencil_end[3];
	
	double t,dt;
	float h, factor, pos[3],x[3],s[3],v[3];
	
	ComputeRHS(const ComputeRHS& c)
	{
		memcpy(this, &c,sizeof(ComputeRHS) );
	}


	ComputeRHS(float dt_, float t_=0): dt(dt_), t(t_){
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = nDim==3?-1:0;
		
		stencil_end[0] = +2;
		stencil_end[1] = +2;
		stencil_end[2] = nDim==3?+2:+1;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
	{
		typedef BlockType B;
		typedef typename BlockType::ElementType E;
		
		double x[3], v[3], h, pos[3], factor;
		int s[3];
		
		pos[0] = info.origin[0];
		pos[1] = info.origin[1]; 
		pos[2] = info.origin[2]; 
		h = info.h[0];
		factor = -1.0/info.h[0];
		
		for(int iz=0; iz<B::sizeZ; iz++)
			for(int iy=0; iy<B::sizeY; iy++)
				for(int ix=0; ix<B::sizeX; ix++)
				{
					x[0] = pos[0] + ix*h;
					x[1] = pos[1] + iy*h;
					x[2] = pos[2] + iz*h;
					
					_velocity<nDim>(x, t+dt, v);
					
					s[0]= v[0]>0? -1 : 0;
					s[1]= v[1]>0? -1 : 0;
					s[2]= v[2]>0? -1 : 0;
					
					if (nDim==3)
						o(ix,iy,iz).tmp = factor*(
												  (i(ix+s[0]+1, iy, iz).evaluate_rho(dt) - i(ix+s[0], iy, iz).evaluate_rho(dt))*v[0] +
												  (i(ix, iy+s[1]+1, iz).evaluate_rho(dt) - i(ix, iy+s[1], iz).evaluate_rho(dt))*v[1] +
												  (i(ix, iy, iz+s[2]+1).evaluate_rho(dt) - i(ix, iy, iz+s[2]).evaluate_rho(dt))*v[2]);
					else
						o(ix,iy,iz).tmp = factor*(
												  (i(ix+s[0]+1, iy, iz).evaluate_rho(dt) - i(ix+s[0], iy, iz).evaluate_rho(dt))*v[0] +
												  (i(ix, iy+s[1]+1, iz).evaluate_rho(dt) - i(ix, iy+s[1], iz).evaluate_rho(dt))*v[1] );
				}
	}
};

struct UpdateRHS
{
	float dt;
	
	UpdateRHS(const UpdateRHS& c):dt(c.dt){}
	UpdateRHS(float dt_): dt(dt_) {}
	
	template<typename BlockType>
	inline void operator()(const BlockInfo& info, BlockType& b) const
	{
		typedef BlockType B;
		typedef typename BlockType::ElementType E;
		
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		E* ptrE = &(b(0));
		
		for(int iE=0; iE<n; iE++, ptrE++)
		{
			ptrE->drho_dt = ptrE->tmp ;
			ptrE->rho -= dt*ptrE->tmp ;
		}
	}
};

struct Integration
{
	float dt;
	
	Integration(const Integration&c):dt(c.dt){}
	Integration(float dt_): dt(dt_) {}
	
	inline void prepare(BlockInfo& b){}
	
	inline void operator()(P& p) const
	{
		p.rho += p.drho_dt*dt;
		p.drho_dt = 0;
	}
	
	template<typename BlockType>
	inline void operator()(const BlockInfo& info, BlockType& b) const
	{
		typedef BlockType B;
		typedef typename BlockType::ElementType E;
		
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		E* ptrE = &(b(0));
		
		for(int iE=0; iE<n; iE++, ptrE++)
		{
			ptrE->rho += ptrE->drho_dt*dt;
			ptrE->drho_dt = 0;
		}
	}
};
	

const float toleranceIC = 0.0001;//0.01;

class Smoke
{
public:
		static const bool bUseSpaceTimeSorter = true;
		static const bool b3D = false;
		static const int nDim = b3D?3:2;
		static const int blocksPerDimension = 4;
		static const int blockSize = 16;
		static const int maxLevel = 5;
		static const int resJump = 1;
		
		typedef Block< P, blockSize, blockSize, b3D? blockSize: 1> B;
		typedef Wavelets_AverageInterp5thOrder W;
		//typedef Wavelets_Interp4thOrder W;
		//typedef Wavelets_AverageInterp3rdOrder W;
#ifdef _MRAG_TBB
		typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;
#else
		typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
#endif
	
		Grid<W, B> grid;
private:
		Refiner_SpaceExtension refiner;
		Compressor compressor;
		BlockLab<B> lab;
		BlockFWT<W, B, smoke_projector> blockfwt;
		GridViewer viewer;
		SpaceTimeSorter stSorter;
		Profiler profiler;
		float t;
		string m_sFormat;
		int m_iCurrentDumpTime;
	
		void _step(float t, float dt)
		{
			BoundaryInfo* boundaryInfo = &grid.getBoundaryInfo();
			
			vector<BlockInfo> vInfo;
			
			if (bUseSpaceTimeSorter)
			{
				stSorter.startSession(dt, 2, 0);
				
				while(true)
				{
					SpaceTimeSorter::ETimeInterval type;
					
					double currTime, currDeltaT;
					int level;
					
					const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);

					if (type == SpaceTimeSorter::ETimeInterval_Start)
					{
						ComputeRHS<nDim> computeRHS(currTime, t);
						UpdateRHS updateRHS(currTime);
				
						BlockProcessing::process<BlockLab>(vInfo, grid.getBlockCollection(), *boundaryInfo, computeRHS);
						BlockProcessing::process(vInfo, grid.getBlockCollection(), updateRHS);
					}
					else
					{
						Integration integration(currTime);
						BlockProcessing::process(vInfo, grid.getBlockCollection(), integration);
					}
					
					if (!bContinue) break;
				}
				
				stSorter.endSession();
			}
			else
			{
				abort();
			}
		}

		template<int nDim>
		static float _ic_func(float x[3])
		{
		
			const float dr = 0.01;
			const float r = sqrt(pow(x[0]-0.35,2) + pow(x[1]-0.35, 2) + (nDim>2?pow(x[2]-0.35,2) : 0));
			
			return min(1., max(0., 1+ (0.15 - r)/dr));
		}
		
		template<int nDim>
		static void _ic(Grid<W, B>& grid)
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
					float x[3];
					
					info.pos(x, ix, iy, iz);
					
					block(ix,iy,iz).rho = _ic_func<nDim>(x);
				}
			}
		}
		
	public:
		Smoke(string sFormat = "Grid_At_Time%03d.%s"):
		grid(blocksPerDimension,blocksPerDimension, b3D?blocksPerDimension:1), 
		refiner(resJump), compressor(resJump), profiler(),
		lab(), blockfwt(),  t(0), viewer(true, true),
		stSorter(), m_sFormat(sFormat), m_iCurrentDumpTime(0)
		{
			
			Profiler prof2;
			prof2.getAgent("Smoke Constructor").start();
			printf("********************* Start Smoke Constructor\n");
			grid.setCompressor(&compressor);
			grid.setRefiner(&refiner);
			
			stSorter.connect(grid);
			

			_ic<nDim>(grid);
			Science::AutomaticRefinement< 0,0 >(grid, blockfwt, toleranceIC/2, maxLevel, -1, &profiler, _ic<nDim>);
			
			_ic<nDim>(grid);
			Science::AutomaticCompression< 0,0 >(grid, blockfwt, toleranceIC/2, -1,   &profiler, _ic<nDim>);
			_ic<nDim>(grid);

			printf("********************* End Smoke Constructor\n");
			prof2.getAgent("Smoke Constructor").stop();
			prof2.printSummary();
			profiler.printSummary();
			profiler.clear();
		}
		
		void Step()
		{
			static int iCounter = 0;
			iCounter++;
			//if (iCounter>4) exit(0);
			
			//for(int k=0; k<10; k++)
			{
				const int nStep = 5;
				const double dx = 1.0/blockSize;
				const double CFL = 0.5;
				const float dt = dx/CFL;//8e-1;
				
				Science::AutomaticRefinement< 0,0 >(grid, blockfwt, toleranceIC/2,maxLevel, -1, &profiler);
				
				{
				
					printf("/////////// MEM MB: %f  (%d blocks) ////////////////// \n", grid.getMemorySize(), grid.getBlocksInfo().size());
				}
				
				printf("Start Computation \n");
				profiler.getAgent("computation").start();
				for(int i=0; i<nStep; i++,t+=dt)
					_step(t, dt);
				profiler.getAgent("computation").stop();
				
				Science::AutomaticCompression< 0,0 >(grid, blockfwt, toleranceIC, -1,  &profiler);
				
				profiler.printSummary();
				//int i; cin>>i;
				//exit(0);
			}
		}
		
		float getTime() { return t;}
		
		void DumpData()
		{
			//1. create suitable filenames
			//2. dump the information of the blocks
			//3. dump the information about the data
			
			//1.
			profiler.getAgent("dumping").start();
			
			char buf[300], buf2[300];
			sprintf(buf, m_sFormat.data(), m_iCurrentDumpTime, "txt");
			sprintf(buf2, m_sFormat.data(), m_iCurrentDumpTime, "grid");
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			//2.
			{
				FILE * file = fopen(buf, "w");
				assert(file!=NULL);
				fprintf (file, "Wavelets: %s", "Wavelets_Interp2ndOrder\n"); 
				fprintf(file, "Cell-centered? %s\n", W::bIsCellCentered ? "yes" : "no");
				fprintf(file, "Blocks: %d\n", vInfo.size());
				fprintf(file, "Block size: %d %d %d\n", blockSize, blockSize, blockSize);
				fprintf(file, "Ghosts at 0: [%d, %d[\n", W::HsSupport[0], W::HsSupport[1]);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					fprintf(file, "Block %d: Tree Index: %d %d %d, %d\n", i,
							info.index[0], info.index[1], info.index[2], info.level); 
				}
				
				fclose(file);
			}
			
			//3.
			{
				FILE * file = fopen(buf2, "wb");
				assert(file!=NULL);
			
				const int steStart[3] ={ -1,-1,-1};
				const int steEnd[3] ={ +2,+2,+2};
				
				lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
				
				const int sX = W::HsSupport[0];
				const int sY = W::HsSupport[0];
				const int sZ = W::HsSupport[0];
				
				const int eX = blockSize + W::HsSupport[1] - 1;
				const int eY = blockSize + W::HsSupport[1] - 1;
				const int eZ = blockSize + W::HsSupport[1] - 1;
				
				Matrix3D<float> * matData = new Matrix3D<float>(eX - sX, eY - sY, eZ - sZ);
				
				for(int i=0; i<vInfo.size(); i++)
				{
					BlockInfo& info = vInfo[i];
					
					lab.load(info);
					
					for(int iz=sZ; iz<eZ; iz++)
					for(int iy=sY; iy<eY; iy++)
					for(int ix=sX; ix<eX; ix++)
						matData->Access(ix-sX, iy-sY, iz-sZ) = (float)lab(ix, iy, iz).rho;
					
					matData->Serialize(file);
				}
				
				fclose(file);
			}
			
			m_iCurrentDumpTime ++;
			profiler.getAgent("dumping").stop();
			
		}
		
		void Render()
		{
			if (!b3D)
				viewer.drawContent(grid, grid.getBlockCollection());

			viewer.drawSketch(grid, false);
			
			if (b3D)
				DumpData();
		}
	};
