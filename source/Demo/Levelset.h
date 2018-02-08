/*
 *  Levelset.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/5/08.
 *  Remodeled by Michael Bergdorf on 9/9/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "../MRAGCore/MRAGCommon.h"
#include "../MRAGCore/MRAGWavelets_Interp2ndOrder.h"
#include "../MRAGCore/MRAGWavelets_Interp4thOrder.h"
#include "../MRAGCore/MRAGBlock.h"
#include "../MRAGCore/MRAGrid.h"
#include "../MRAGCore/MRAGBlockFWT.h"
#include "../MRAGScience/MRAGScienceCore.h"
#include "../MRAGScience/MRAGSpaceTimeSorter.h"
#include "../MRAGscience/MRAGSimpleLevelsetBlock.h"
#include "../MRAGscience/MRAGRefiner_SpaceExtension.h"


using namespace MRAG;

struct PLS
{
	float rho;
	float drho_dt;
   float tmp;
	PLS(): rho(0), drho_dt(0), tmp(0) {}
	
	PLS(float rho_, float drho_dt_): rho(rho_), drho_dt(drho_dt_) {}
	
   Real levelset() const{
      return rho;
   }
   
	void operator += (PLS t)
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
      tmp = 0;
	}
	
	float evaluate_rho(double dt)
	{
		return rho + drho_dt*dt;
	}
};

PLS operator*(const PLS& p, Real v)
{
	PLS t(p.rho*v, p.drho_dt*v);
	return t;
}

/*
template <typename T, int i> inline Real levelset_projector_impl(const T&t)
{
	return (Real)t.rho;
}
*/

make_projector(levelset_projector, SimpleLevelsetBlock<PLS>::levelset_projector_impl)



const float toleranceIC = 0.0001;

class Levelset
	{
		static const bool bUseSpaceTimeSorter = true;
		static const bool b3D = false;
		static const int blocksPerDimension = 1;
		static const int blockSize = 8;
		static const int maxLevel = 6;
		static const int resJump = 2;
		static const int narrowBandWidth = 6;
		
		typedef SimpleLevelsetBlock< PLS, narrowBandWidth, blockSize, blockSize, b3D? blockSize: 1> B;
		typedef Wavelets_Interp2ndOrder W;
		
		Grid<W, B> grid;
		Refiner_SpaceExtension refiner;
		Compressor compressor;
		BlockLab<B> lab;
		BlockFWT<W, B, levelset_projector> blockfwt;
		GridViewer viewer;
		SpaceTimeSorter stSorter;
		float t;
		
		void _step(float t, float dt)
		{
			vector<BlockInfo> vInfo;
			if (bUseSpaceTimeSorter)
			{
				double currTime, currDeltaT;
				int level;
				stSorter.startSession(dt, 2, 0);
				
				while(true)
				{
					SpaceTimeSorter::ETimeInterval type;
					const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
					//printf("level =%d, t = %f, dt= %f\n", level , currTime, currDeltaT);
					
					if (type == SpaceTimeSorter::ETimeInterval_Start)
					{
						//	printf("COMPUTE RHS %f\n", currDeltaT);
						//					double time, double currentDT, double dt
						_step_computeRHS_TR(vInfo, currTime );
					}
					else
					{
						//printf("INTEGRATE %f\n", currDeltaT);
						_step_integrate(vInfo, currTime);
					}
					
					if (!bContinue) break;
				}
				
				stSorter.endSession();
			}
			else
			{
				abort();
				vInfo = grid.getBlocksInfo();
				
				_step_computeRHS(vInfo);
				_step_integrate(vInfo,dt);
			}
		}
		

		void _step_computeRHS_TR(vector<BlockInfo>& vInfo,  double dt)
		{
			const int steStart[3] ={ -1,-1,0};
			const int steEnd[3] ={ +2,+2,+1};
			
			lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				lab.load(info);
				
				
				const float pos[3] = {info.origin[0], info.origin[1], info.origin[3]};
				const float h = info.h[0];
				const float factor = -1.0/info.h[0];
				
				float v[3], x[3];
				int s[2];
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
					{
						x[0] = pos[0] + ix*h;
						x[1] = pos[1] + iy*h;
						x[2] = 0;
						
						_velocity(x, t, v);
						
						s[0]= v[0]>0? -1 : 0;
						s[1]= v[1]>0? -1 : 0;
						
						block(ix,iy).tmp = factor*(
                                             (lab(ix+s[0]+1, iy).evaluate_rho(dt) - lab(ix+s[0], iy).evaluate_rho(dt))*v[0] +
                                             (lab(ix, iy+s[1]+1).evaluate_rho(dt) - lab(ix, iy+s[1]).evaluate_rho(dt))*v[1]);
					}
			}
      	
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++){
						block(ix,iy).drho_dt = block(ix,iy).tmp ;
                  block(ix,iy).rho -= dt*block(ix,iy).tmp;                  
               }
			}		
		}
		
         
		float _ic_func(float x[3])
		{
			const float r = sqrt(pow(x[0]-0.5,2) + pow(x[1]-0.70, 2));
			return -r+0.15;
		}
		
		void _ic()
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				block.setH(info.h[0]);
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							float x[3];
							
							info.pos(x, ix, iy, iz);
							
							block(ix,iy,iz).rho = _ic_func(x);
						}
			}
		}
		
		void _velocity(float x[3], float t, float v[3])
		{
			const double p[3] = { M_PI*(x[0] - 0.00), M_PI*(x[1] - 0.00), M_PI*(x[2] - 0.00)};
			const double factor = 0.1;
			v[0] = -factor*2.0* pow(sin(p[0]), 2) * sin(p[1]) * cos(p[1]);
			v[1] = factor*2.0* pow(sin(p[1]), 2) * sin(p[0]) * cos(p[0]);
			v[2] = 0;
		}
		
		
		
		
		void _step_computeRHS(vector<BlockInfo>& vInfo)
		{
			const int steStart[3] ={ -1,-1,0};
			const int steEnd[3] ={ +2,+2,+1};
			
			lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				lab.load(info);
				
				
				const float pos[3] = {info.origin[0], info.origin[1], info.origin[3]};
				const float h = info.h[0];
				const float factor = -1.0/info.h[0];
				
				float v[3], x[3];
				int s[2];
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
					{
						x[0] = pos[0] + ix*h;
						x[1] = pos[1] + iy*h;
						x[2] = 0;
						
						_velocity(x, t, v);
						
						s[0]= v[0]>0? -1 : 0;
						s[1]= v[1]>0? -1 : 0;
						
						block(ix,iy).drho_dt = factor*(
                                                 (lab(ix+s[0]+1, iy).rho - lab(ix+s[0], iy).rho)*v[0] +
                                                 (lab(ix, iy+s[1]+1).rho - lab(ix, iy+s[1]).rho)*v[1]);
					}
			}
		}
		
		void _step_integrate(vector<BlockInfo>& vInfo, float dt)
		{
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
						block(ix,iy).integrate(dt);
			}		
		}
		
		/*void _step_bfecc(float t, float dt)
		 {
		 printf("BFECC!\n");
		 vector<BlockInfo> vInfo = grid.getBlocksInfo();
		 
		 _step_computeRHS(vInfo);
		 
		 for(int i=0; i<vInfo.size(); i++)
		 {
		 BlockInfo& info = vInfo[i];
		 B& block = grid.getBlockCollection()[info.blockID];
		 
		 for(int iy=0; iy<B::sizeY; iy++)
		 for(int ix=0; ix<B::sizeX; ix++)
		 {
		 block(ix,iy).rho_tilda = block(ix,iy).rho;
		 block(ix,iy).rho += dt*block(ix,iy).drho_dt;
		 }
		 }
		 
		 _step_computeRHS(vInfo);
		 
		 for(int i=0; i<vInfo.size(); i++)
		 {
		 BlockInfo& info = vInfo[i];
		 B& block = grid.getBlockCollection()[info.blockID];
		 float tmp;
		 for(int iy=0; iy<B::sizeY; iy++)
		 for(int ix=0; ix<B::sizeX; ix++)
		 {
		 tmp = block(ix,iy).rho_tilda;
		 
		 block(ix,iy).rho_tilda = block(ix,iy).rho;
		 block(ix,iy).rho_bar = block(ix,iy).rho -dt*block(ix,iy).drho_dt;
		 block(ix,iy).rho = tmp*1.5 - 0.5*block(ix,iy).rho_bar;
		 }
		 }
		 
		 
		 _step_computeRHS(vInfo);
		 
		 _step_integrate(vInfo,t,dt);
		 
		 }*/
      

      
      void _drawGridPoint3f(const int block_index[3],  int level, const float point_index[3], const float *vColorBoundary = NULL, double factor=1.)
      {
         const double h[2]= {pow(2.,-level), pow(2.,-level)};
         const bool bCollocated = true;
         const double start[2] = {block_index[0]*h[0],block_index[1]*h[1]};
         //const double end[2] = {(block_index[0]+1)*h[0],(block_index[1]+1)*h[1]};
         const int n[2] = {B::sizeX+1, B::sizeY+1};
         const double d[2] = {h[0]/(n[0]-1),h[1]/(n[1]-1)};
         const double point_center[2] = {start[0] + (point_index[0]+ (bCollocated?0: 0.5))*d[0], start[1] + (point_index[1]+(bCollocated?0: 0.5))*d[1]};
         
         glPointSize(600*d[0]/2*factor);
         glBegin(GL_POINTS);
         if (vColorBoundary != NULL)
            glColor3fv(vColorBoundary);
         else
            glColor3f(1,0,0);
         
         
         glVertex2f(point_center[0], point_center[1]);
         glEnd();
      }
      
      void _drawLevelsetIntersections(){
         const float vContent[3] = {1.0,0.0,0.573};

         vector<BlockInfo> vInfo = grid.getBlocksInfo();
         for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
            const int level = info.level;
            const int block_index[] = {
               info.index[0],info.index[1],info.index[2]
            };
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY-1; iy++)
						for(int ix=0; ix<B::sizeX-1; ix++)
						{
                     if(block(ix,iy,iz).rho*block(ix+1,iy,iz).rho<0.0){
                        const float x[] = {
                           (float)ix + block(ix,iy,iz)/(block(ix,iy,iz)-block(ix+1,iy,iz)),
                           (float)iy,
                           0.0
                        };
                        _drawGridPoint3f(block_index, level, x,vContent,1.0);
                     } 
                     if(block(ix,iy,iz).rho*block(ix,iy+1,iz).rho<0.0){
                        const float x[] = {
                           (float)ix,
                           (float)iy + block(ix,iy,iz)/(block(ix,iy,iz)-block(ix,iy+1,iz)),
                           0.0
                        };
                        _drawGridPoint3f(block_index, level, x,vContent,1.0);
                     } 
						}
         }
      }
	public:
		Levelset():
		grid(blocksPerDimension,blocksPerDimension, b3D?blocksPerDimension:1), 
		refiner(resJump), compressor(resJump), 
		lab(), blockfwt(), viewer(), t(0),
		stSorter()
		{
			grid.setCompressor(&compressor);
			grid.setRefiner(&refiner);
			
			stSorter.connect(grid);
			
			_ic();
         Science::AutomaticRefinementForLevelsets(grid, blockfwt, toleranceIC, maxLevel);
			_ic();
         Science::AutomaticCompressionForLevelsets(grid, blockfwt, toleranceIC);         
			_ic();
		}
		
		void Step()
		{
			//for(int k=0; k<10; k++)
			{
				const int nStep = 1;
				const double dx = 1.0/blockSize;
				const double CFL = 0.2;
				const float dt = dx/CFL;//8e-1;
				
            Science::AutomaticRefinementForLevelsets(grid, blockfwt, toleranceIC/2, maxLevel);

				for(int i=0; i<nStep; i++)
					_step(t+=dt, dt);

            Science::AutomaticCompressionForLevelsets(grid, blockfwt, toleranceIC);         

			}
		}
		
		float getTime() { return t;}
		

		
		void Render()
		{
			viewer.drawContent(grid, grid.getBlockCollection());
			viewer.drawSketch(grid,false);
         _drawLevelsetIntersections();

         
		}
	};
