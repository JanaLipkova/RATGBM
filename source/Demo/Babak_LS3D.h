/*
 *  Babak_LS3D.h
 *  MRAG
 *
 *  Created by Babak Hejazialhosseini  on 11/5/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


/*#include "../MRAGCore/MRAGCommon.h"
 #include "../MRAGCore/MRAGWavelets_Interp2ndOrder.h"
 #include "../MRAGCore/MRAGWavelets_Interp4thOrder.h"
 #include "../MRAGCore/MRAGBlock.h"
 #include "../MRAGCore/MRAGrid.h"
 #include "../MRAGCore/MRAGBlockFWT.h"
 #include "../MRAGScience/MRAGScienceCore.h"
 #include "../MRAGScience/MRAGSpaceTimeSorter.h"
 #include "../MRAGscience/MRAGSimpleLevelsetBlock.h"
 #include "../MRAGscience/MRAGRefiner_SpaceExtension.h"
 #include "../MRAGvisual/GridViewer.h"
 #include "../MRAGvisual/MRAGVisualTypes.h"*/

#include "../MRAGcore/MRAGCommon.h"
#include "../MRAGcore/MRAGWavelets_Interp2ndOrder.h"
#include "../MRAGcore/MRAGWavelets_Interp4thOrder.h"
#include "../MRAGcore/MRAGBlock.h"
#include "../MRAGcore/MRAGrid.h"
#include "../MRAGcore/MRAGBlockFWT.h"
#include "../MRAGscience/MRAGScienceCore.h"
#include "../MRAGscience/MRAGSpaceTimeSorter.h"
#include "../MRAGscience/MRAGRefiner_SpaceExtension.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_SingleCPU.h"
#include "../MRAGmultithreading/MRAGBlockProcessing_TBB.h"
#include "../MRAGcore/MRAGProfiler.h"
#include "../MRAGvisual/MRAGVisualTypes.h"
#include "../MRAGvisual/GridViewer.h"
#include "../MRAGio/MRAG_IO_Native.h"
#include "../MRAGio/MRAG_IO_VTK.h"

using namespace MRAG;


struct PLS
{
	float rho,rho_0;
	float drho_dt;
	float tmp;
	PLS(): rho(0), drho_dt(0), tmp(0) {}
	
	PLS(float rho_, float drho_dt_): rho(rho_), drho_dt(drho_dt_) {}
	PLS(float rho_): rho(rho_), drho_dt(0) {}

	Real levelset() const{
		return rho;
	}
	
	void operator += (PLS t) 
	{
		rho_0 +=t.rho_0;
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


 template <typename T, int i> inline Real levelset_projector_impl(const T&t)
 {
 return (Real)t.rho;
 }
 

make_projector(levelset_projector, SimpleLevelsetBlock<PLS>::levelset_projector_impl)



const float toleranceIC = 0.001;

const int maxStencil[2][3] = {
-3,-3,-3,
4,4,4
};

class Levelset
	{
	public:
		static const bool bRestartFromFile = true;
		static const bool bPostProcessData = false;
		static const bool bUseSpaceTimeSorter = true;
		static const int nStepsPerSerialization = 4;
		static const bool b3D = true;
		static const int blocksPerDimension = 1;
		static const int blockSize = 26;
		static const int maxLevel = 4;
		static const int resJump = 2;
		static const int narrowBandWidth = 6;
		static const int nDim = b3D?3:2;
		
		typedef SimpleLevelsetBlock< PLS, narrowBandWidth, blockSize, blockSize, b3D? blockSize: 1> B;
		//typedef Wavelets_Interp2ndOrder W;
		typedef Wavelets_Interp4thOrder W;
		Grid<W, B> grid;
		

	private:	
		Refiner_SpaceExtension refiner;
		Compressor compressor;
		BlockLab<B> lab;
		BlockFWT<W, B, levelset_projector> blockfwt;
		GridViewer viewer;
		SpaceTimeSorter stSorter;
		Profiler profiler;
		string m_sFormat, m_sToRenderFormat;
		int m_iCurrentDumpTime;
		float t;
		IO_VTK< W, B, levelset_projector > vtk;
		
		//typedef Multithreading::BlockProcessing_SingleCPU<B> BlockProcessing;
		typedef Multithreading::BlockProcessing_TBB<B> BlockProcessing;
		
		void _step(float t, float dt)
		{
			char fname[20];
			
			if (bUseSpaceTimeSorter)
			{
				
				BoundaryInfo* boundaryInfo = &grid.getBoundaryInfo();
								
				stSorter.startSession(dt, 2, 0);
					
				if(t>9.0) exit(0);				
				
				while(true)
				{
					double currTime, currDeltaT;
					int level;
					SpaceTimeSorter::ETimeInterval type;
					vector<BlockInfo> vInfo;

					const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
					
					if (type == SpaceTimeSorter::ETimeInterval_Start)
					{
						//BP_ComputeRHS_TR task1(currTime);

						BP_ComputeRHS_TR_WENO5 task1(currTime, (t+currTime<4.5)? 1 : -1);
						BP_UpdateRHS_TR task2(currTime);

						BlockProcessing::process< BlockLab >(vInfo, grid.getBlockCollection(), *boundaryInfo, task1);
						BlockProcessing::process(vInfo, grid.getBlockCollection(), task2);
					}
					else
					{
						BP_Integrate_TR task(currTime);
						BlockProcessing::process(vInfo, grid.getBlockCollection(), task);
						
					}
					
					if (!bContinue) break;
				}
				
				//sprintf(fname,"%s%i","babak",int(100*t));
				//vtk.Write(grid, fname);
				cout<<"I am at time= "<<t<<endl;
				printf("/////////// MEM MB: %f  (%d blocks) ////////////////// \n", grid.getMemorySize(), grid.getBlocksInfo().size());
				
				// reinitialize here
				reinit_run();
				
				stSorter.endSession();
			}
			else
			{
				abort();
			}
		}
		
		void RefineOnceEverywhere(int maxResJump)
		{
			for(int i=0; i<=maxResJump; i++)
			{
				blockfwt.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
				vector<BlockInfo> vInfo = grid.getBlocksInfo();
				
				set<int> shouldBeRefined;
				for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
				{
					blockfwt.fwt<0>(*it);
					
					if (blockfwt.getReport().getOverAll_DetailMaxMag()>0)
						shouldBeRefined.insert(it->blockID);
				}
				
				grid.refine(shouldBeRefined);
			}
			
			
			/*	{
			 vector<BlockInfo> vInfo = grid.getBlocksInfo();
			 blockfwt.prepare(grid.getBlockCollection(), grid.getBoundaryInfo());
			 set<int> shouldBeRefined;
			 for(vector<BlockInfo>::const_iterator it = vInfo.begin(); it != vInfo.end(); it++)
			 {
			 blockfwt.fwt<0>(*it);
			 assert(blockfwt.getReport().getOverAll_DetailMaxMag() < 1e-15);
			 }
			 }*/
			
		}
		
		template<typename RealType>
		static void _velocity(const RealType x[3],  RealType  t, RealType  v[3])
		{
			//const double p[3] = { M_PI*(x[0] - 0.00), M_PI*(x[1] - 0.00), M_PI*(x[2] - 0.00)};
			const double factor = 0.2;
			
			v[0] =  factor*2.0* pow(sin(M_PI*x[0]), 2) * sin(2*M_PI*x[1]) * sin(2*M_PI*x[2]);
			v[1] = -factor* pow(sin(M_PI*x[1]), 2) * sin(2*M_PI*x[0]) *  sin(2*M_PI*x[2]);
			v[2] = -factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1]);
			
			
		}
		
		template<typename RealType>
		static void _velocity2(const RealType x[3], RealType  t,RealType  v[3])
		{
			//const double p[3] = { M_PI*(x[0] - 0.00), M_PI*(x[1] - 0.00), M_PI*(x[2] - 0.00)};
			const double factor = 0.2;
			
			v[0] =  factor*2.0* pow(sin(M_PI*x[0]), 2) * sin(2*M_PI*x[1]) * sin(2*M_PI*x[2]);
			v[1] = -factor* pow(sin(M_PI*x[1]), 2) * sin(2*M_PI*x[0]) *  sin(2*M_PI*x[2]);
			v[2] = -factor* pow(sin(M_PI*x[2]), 2) * sin(2*M_PI*x[0]) * sin(2*M_PI*x[1]);
			
			
			v[0]=-1.0*v[0]; 
			v[1]=-1.0*v[1];
			v[2]=-1.0*v[2];
			
			
		}
		
		struct BP_UpdateRHS_TR
		{
			double dt;

			BP_UpdateRHS_TR(const double dt_): dt(dt_){}
			BP_UpdateRHS_TR(const BP_UpdateRHS_TR& c): dt(c.dt) {}

			template<typename B>
			inline void operator()(BlockInfo& info, B& b) const
			{
				typedef typename B::ElementType E;

				const int n = B::sizeZ*B::sizeY*B::sizeX;
		
				E* ptrE = &(b[0]);

				for(int iE=0; iE<n; iE++, ptrE++)
				{
					ptrE->drho_dt = ptrE->tmp;
					ptrE->rho -= dt*ptrE->tmp;
				}		
			}
		};

				
		struct BP_Integrate_TR
		{
			double dt;

			BP_Integrate_TR(const double dt_): dt(dt_){}
			BP_Integrate_TR(const BP_Integrate_TR& c): dt(c.dt) {}

			template<typename B>
			inline void operator()(BlockInfo& info, B& b) const
			{
				typedef typename B::ElementType E;

				const int n = B::sizeZ*B::sizeY*B::sizeX;
		
				E* ptrE = &(b[0]);

				for(int iE=0; iE<n; iE++, ptrE++)
				{
					ptrE->rho += dt*ptrE->drho_dt;
					ptrE->drho_dt = 0;
				}		
			}
		};

	
		struct BP_ComputeRHS_TR
		{
			int stencil_start[3], stencil_end[3];
			Real t,dt;
			
			BP_ComputeRHS_TR(const BP_ComputeRHS_TR& c)
			{
				memcpy(this, &c, sizeof(BP_ComputeRHS_TR) );
			}
			
			BP_ComputeRHS_TR(float dt_, float t_=0): dt(dt_), t(t_){
				stencil_start[0] = -1;
				stencil_start[1] = -1;
				stencil_start[2] = -1;
				
				stencil_end[0] = +2;
				stencil_end[1] = +2;
				stencil_end[2] = +2;
			}

			inline const Real Upwind_Flux(
				const Real u_plus, const Real u_minus, 
				const Real v_plus, const Real v_minus, 
				const Real w_plus, const Real w_minus, const float v[3]) const
			{
				
				Real term1=0,term2=0,term3=0;
							
				if (v[0]<0)
					term1 = u_plus;
				else if (v[0]>0)
					term1 = u_minus;
				else
					term1 = 0.5*(u_minus+u_plus);
				
				
				if (v[1]<0)
					term2 = v_plus;
				else if (v[1]>0)
					term2 = v_minus;
				else
					term2 = 0.5*(v_minus+v_plus);
				
				
				if (v[2]<0)
					term3 = w_plus;
				else if (v[2]>0)
					term3 = w_minus;
				else
					term3 = 0.5*(w_minus+w_plus);
				
				return (v[0]*term1 + v[1]*term2 + v[2]*term3);
			}
			
			template <typename Field>
			const Real mainRHS(Field& f, const Real spacing, const int ix, const int iy, const int iz, const float v[3], const float dt) const
			{
				const Real u_minus = (-f(ix-1,iy, iz).evaluate_rho(dt) + f(ix, iy,iz).evaluate_rho(dt))/spacing;
				const Real u_plus  = (-f(ix,iy, iz).evaluate_rho(dt) + f(ix+1, iy,iz).evaluate_rho(dt))/spacing;
				const Real v_minus = (-f(ix,iy-1, iz).evaluate_rho(dt) + f(ix, iy,iz).evaluate_rho(dt))/spacing;
				const Real v_plus  = (-f(ix,iy, iz).evaluate_rho(dt) + f(ix, iy+1,iz).evaluate_rho(dt))/spacing;
				const Real w_minus = (-f(ix,iy, iz-1).evaluate_rho(dt) + f(ix, iy,iz).evaluate_rho(dt))/spacing;
				const Real w_plus  = (-f(ix,iy, iz).evaluate_rho(dt) + f(ix, iy,iz+1).evaluate_rho(dt))/spacing;
				
				return Upwind_Flux( u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, v);
			}

			template<typename LabType, typename BlockType>
			inline void operator()(LabType& i, BlockInfo& info, BlockType& o) const
			{
				typedef BlockType B;
				typedef typename BlockType::ElementType E;
				
				float x[3], v[3], t=0;
				int s[3];
				
				const Real pos[3] = {
					info.origin[0],
					info.origin[1], 
					info.origin[2]
				}; 

				const Real h = info.h[0];

				const Real factor = -1.0/info.h[0];
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							x[0] = pos[0] + ix*h;
							x[1] = pos[1] + iy*h;
							x[2] = pos[2] + iz*h;
							
							Levelset::_velocity(x, (float)(t+dt), v);

							o(ix,iy,iz).tmp = -mainRHS(i, h, ix, iy ,iz, v, dt);
						}		
			}
		};
				
		struct BP_ComputeRHS_TR_WENO5
		{
			int stencil_start[3], stencil_end[3];
			Real dt, sign_val;

			BP_ComputeRHS_TR_WENO5(const BP_ComputeRHS_TR_WENO5& c)
			{
				memcpy(this, &c, sizeof(BP_ComputeRHS_TR_WENO5) );
			}

			BP_ComputeRHS_TR_WENO5(Real dt_, Real sign_): dt(dt_), sign_val(sign_)
			{
				stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
				stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
			}

			template<typename LabType, typename BlockType>
			inline void operator()(LabType& i, BlockInfo& info, BlockType& o) const
			{	
				typedef BlockType B;
				typedef typename BlockType::ElementType E;
				
				const Real h = info.h[0];
								
				Real x[3], v[3];
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{		
							info.pos(x,ix,iy,iz);
							_velocity(x, 0.0f, v);

							v[0] *= sign_val;
							v[1] *= sign_val;
							v[2] *= sign_val;

							o(ix,iy,iz).tmp = -mainRHS(i, h , ix, iy, iz, v, dt);
						}
			}
			
			inline const Real Upwind_Flux(const Real u_plus, const Real u_minus, const Real v_plus, const Real v_minus, const Real w_plus, const Real w_minus, const Real v[3]) const
			{
				const Real term[3] = {
					(v[0]<0)? u_plus : ((v[0]>0)? u_minus : 0.5*(u_minus+u_plus)),
					(v[1]<0)? v_plus : ((v[1]>0)? v_minus : 0.5*(v_minus+v_plus)),
					(v[2]<0)? w_plus : ((v[2]>0)? w_minus : 0.5*(w_minus+w_plus)),
				};
				
				return v[0]*term[0] + v[1]*term[1] + v[2]*term[2];
			}
			
			inline const Real phi_WENO(const Real a, const Real b, const Real c, const Real d) const
			{
				const Real eps = 1e-6;
				
				const Real IS[3] = {
					13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b),
					13*(b-c)*(b-c) + 3*(b+c)*(b+c),
					13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d)
				};
				
				const Real alpha[3] = {
					1./((eps + IS[0])*(eps + IS[0])),
					6./((eps + IS[1])*(eps + IS[1])),
					3./((eps + IS[2])*(eps + IS[2])),
				};
				
				const Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
				
				const Real w[2] = {
					alpha[0]/sum_alpha,
					alpha[2]/sum_alpha
				};
				
				return 1./3*w[0]*(a-2*b+c) + 1./6*(w[1]-0.5)*(b-2*c+d);
			}
			
			template <int index, typename Field> inline const Real Dplus(Field& f, const int ix, const int iy, const int iz) const
			{ 
				//return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).rho - f(ix, iy, iz).rho);
				return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).evaluate_rho(dt) - f(ix, iy, iz).evaluate_rho(dt));
			}
			
			template <int index, typename Field> inline const Real Dminus(Field& f, const int ix, const int iy, const int iz) const
			{ 
				//return (f(ix, iy, iz).rho - f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).rho );
				return (f(ix, iy, iz).evaluate_rho(dt) - f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).evaluate_rho(dt) );
			}
			
			template <int index, typename Field> inline const Real DplusDminus(Field& f, const int ix, const int iy, const int iz) const
			{ 
				return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).evaluate_rho(dt) +
						f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).evaluate_rho(dt) +
						-2*f(ix, iy,iz).evaluate_rho(dt));
			}
			
			template <typename Field>
			inline const Real mainRHS(Field& f, const Real spacing, const int ix, const int iy, const int iz, const Real v[3], const Real dt) const
			{
				const Real u_common_term = 1./(12*spacing)*
				(-Dplus<0>(f, ix-2, iy, iz) + 7*Dplus<0>(f,ix-1, iy, iz) + 7*Dplus<0>(f,ix, iy, iz) -Dplus<0>(f, ix+1, iy, iz) );
				
				const Real v_common_term = 1./(12*spacing)*
				(-Dplus<1>(f, ix, iy-2, iz) + 7*Dplus<1>(f,ix, iy-1, iz) + 7*Dplus<1>(f,ix, iy, iz) -Dplus<1>(f, ix, iy+1, iz) );
				
				const Real w_common_term = 1./(12*spacing)*
				(-Dplus<2>(f, ix, iy, iz-2) + 7*Dplus<2>(f,ix, iy, iz-1) + 7*Dplus<2>(f,ix, iy, iz) -Dplus<2>(f, ix, iy, iz+1) );
				
				const Real phi_weno_x[2] = {
					phi_WENO(DplusDminus<0>(f, ix-2, iy, iz)/spacing,	DplusDminus<0>(f, ix-1, iy, iz)/spacing, 
							 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix+1, iy, iz)/spacing),
					phi_WENO(DplusDminus<0>(f, ix+2, iy, iz)/spacing,	DplusDminus<0>(f, ix+1, iy, iz)/spacing, 
							 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix-1, iy, iz)/spacing)
					
				};
				
				const Real phi_weno_y[2] = {
					phi_WENO(DplusDminus<1>(f, ix, iy-2, iz)/spacing,	DplusDminus<1>(f, ix, iy-1, iz)/spacing, 
							 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy+1, iz)/spacing),
					phi_WENO(DplusDminus<1>(f, ix, iy+2, iz)/spacing,	DplusDminus<1>(f, ix, iy+1, iz)/spacing, 
							 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy-1, iz)/spacing)
				};
				
				const Real phi_weno_z[2] = {
					phi_WENO(DplusDminus<2>(f, ix, iy, iz-2)/spacing,	DplusDminus<2>(f, ix, iy, iz-1)/spacing, 
							 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz+1)/spacing),
					phi_WENO(DplusDminus<2>(f, ix, iy, iz+2)/spacing,	DplusDminus<2>(f, ix, iy, iz+1)/spacing, 
							 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz-1)/spacing)
				};
				
				const Real u_minus = u_common_term - phi_weno_x[0];
				const Real u_plus  = u_common_term + phi_weno_x[1];
				const Real v_minus = v_common_term - phi_weno_y[0];
				const Real v_plus  = v_common_term + phi_weno_y[1];
				const Real w_minus = w_common_term - phi_weno_z[0];
				const Real w_plus  = w_common_term + phi_weno_z[1];
				
				return Upwind_Flux( u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, v);
			}
		};
	
		struct BP_ComputeRHS_TR_REV
		{
			int stencil_start[3], stencil_end[3];
			Real dt;
			
			BP_ComputeRHS_TR_REV(const BP_ComputeRHS_TR_REV& c)
			{
				memcpy(this, &c, sizeof(BP_ComputeRHS_TR_REV) );
			}
			
			BP_ComputeRHS_TR_REV(Real dt_): dt(dt_){
				stencil_start[0] = -1;
				stencil_start[1] = -1;
				stencil_start[2] = -1;
				
				stencil_end[0] = +2;
				stencil_end[1] = +2;
				stencil_end[2] = +2;
			}
			
			template<typename LabType, typename BlockType>
			inline void operator()(LabType& i, BlockInfo& info, BlockType& o) const
			{
				typedef BlockType B;
				typedef typename BlockType::ElementType E;
				
				Real x[3], v[3], t=0;
				int s[3];
				
				const Real pos[3] = {
					info.origin[0],
					info.origin[1], 
					info.origin[2]
				}; 

				const Real h = info.h[0];

				const Real factor = -1.0/info.h[0];
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
						{
							x[0] = pos[0] + ix*h;
							x[1] = pos[1] + iy*h;
							x[2] = pos[2] + iz*h;
							
							Levelset::_velocity2(x, t+dt, v);
							
							s[0]= v[0]>0? -1 : 0;
							s[1]= v[1]>0? -1 : 0;
							s[2]= v[2]>0? -1 : 0;
							
							
							o(ix,iy,iz).tmp = factor*(
													  (i(ix+s[0]+1, iy, iz).evaluate_rho(dt) - i(ix+s[0], iy, iz).evaluate_rho(dt))*v[0] +
													  (i(ix, iy+s[1]+1, iz).evaluate_rho(dt) - i(ix, iy+s[1], iz).evaluate_rho(dt))*v[1] +
													  (i(ix, iy, iz+s[2]+1).evaluate_rho(dt) - i(ix, iy, iz+s[2]).evaluate_rho(dt))*v[2]);
						}		
			}
		};
		
		void _n_step_computeRHS_TR(vector<BlockInfo>& vInfo,  double dt)
		{
			const int steStart[3] ={ -1,-1,0};
			const int steEnd[3] ={ +2,+2,+1};
			
			lab.prepare(grid.getBlockCollection(), grid.getBoundaryInfo(),steStart,steEnd);
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				lab.load(info);
				
				
				const float pos[3] = {info.origin[0], info.origin[1], info.origin[2]};
				const float h = info.h[0];
				const float factor = -1.0/info.h[0];
				
				float v[3], gradp[3],gradn[3], sgn, gradmag,igradmag;
				int s[3];
				for(int iz=0; iz<B::sizeZ; iz++)
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
					{
						
						gradp[0]=(lab(ix+1, iy, iz).evaluate_rho(dt) - lab(ix, iy, iz).evaluate_rho(dt))/h;
						gradp[1]=(lab(ix, iy+1, iz).evaluate_rho(dt) - lab(ix, iy, iz).evaluate_rho(dt))/h;
						gradp[2]=(lab(ix, iy, iz+1).evaluate_rho(dt) - lab(ix, iy, iz).evaluate_rho(dt))/h;
						
						gradn[0]=(lab(ix, iy, iz).evaluate_rho(dt) - lab(ix-1, iy, iz).evaluate_rho(dt))/h;
						gradn[1]=(lab(ix, iy, iz).evaluate_rho(dt) - lab(ix, iy-1, iz).evaluate_rho(dt))/h;
						gradn[2]=(lab(ix, iy, iz).evaluate_rho(dt) - lab(ix, iy, iz-1).evaluate_rho(dt))/h;
						
						if (lab(ix, iy, iz).evaluate_rho(dt) > 0.0)
						{ 
							sgn=1.0;
						}
						else
						{
							if (lab(ix, iy, iz).evaluate_rho(dt) < 0.0)
							{
								sgn = -1.0;
							}
							else
							{
								sgn = 0.0;
							}
							
						}
						
						if ((gradp[0]*sgn < 0.0)&&((gradp[0]+gradn[0])*sgn< 0.0))
						{
							v[0] = gradp[0];
						}
						else
						{
							if ((gradn[0]*sgn > 0.0)&&((gradp[0]+gradn[0])*sgn> 0.0)) 
							{
								v[0] = gradn[0];
							}
							else
							{
								v[0] = 0.5* (gradp[0]+gradn[0]);
							}
						}
						
						if ((gradp[1]*sgn < 0.0)&&((gradp[1]+gradn[1])*sgn< 0.0))
						{
							v[1] = gradp[1];
						}
						else
						{
							if ((gradn[1]*sgn > 0.0)&&((gradp[1]+gradn[1])*sgn> 0.0)) 
							{
								v[1] = gradn[1];
							}
							else
							{
								v[1] = 0.5* (gradp[1]+gradn[1]);
							}
						}
						
						if ((gradp[2]*sgn < 0.0)&&((gradp[2]+gradn[2])*sgn< 0.0))
						{
							v[2] = gradp[2];
						}
						else
						{
							if ((gradn[2]*sgn > 0.0)&&((gradp[2]+gradn[2])*sgn> 0.0)) 
							{
								v[2] = gradn[2];
							}
							else
							{
								v[2] = 0.5* (gradp[2]+gradn[2]);
							}
						}
						
						
						gradmag = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
						igradmag = 1.0/gradmag;
						if (gradmag<0.001) igradmag = 0.0;
						
						v[0] = 0.05*v[0] * igradmag;
						v[1] = 0.05*v[1] * igradmag;
						v[2] = 0.05*v[2] * igradmag;
						
						s[0]= v[0]>0? -1 : 0;
						s[1]= v[1]>0? -1 : 0;
						s[2]= v[2]>0? -1 : 0;
						
						block(ix,iy,iz).tmp = factor*(
												   (lab(ix+s[0]+1, iy, iz).evaluate_rho(dt) - lab(ix+s[0], iy, iz).evaluate_rho(dt))*v[0] +
												   (lab(ix, iy+s[1]+1, iz).evaluate_rho(dt) - lab(ix, iy+s[1], iz).evaluate_rho(dt))*v[1] +
												   (lab(ix, iy, iz+s[2]+1).evaluate_rho(dt) - lab(ix, iy, iz+s[2]).evaluate_rho(dt))*v[2]);
					}
			}
			
			for(int i=0; i<vInfo.size(); i++)
			{
				BlockInfo& info = vInfo[i];
				B& block = grid.getBlockCollection()[info.blockID];
				
				for (int iz=0; iz<B::sizeZ; iz++)
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++){
						block(ix,iy,iz).drho_dt = block(ix,iy,iz).tmp ;
						block(ix,iy,iz).rho -= dt*block(ix,iy,iz).tmp;                  
					}
			}		
		}
	
		float _ic_func(float x[3])
		{
			const float r = sqrt(pow(x[0]-0.5,2) + pow(x[1]-0.70, 2) + pow(x[2]-0.70, 2));
			return -r+0.15;
		}
		
		float _ic_func4(float x[3]) const
		{
			const float dr = 0.01;
			const float r = sqrt(pow(x[0]-0.35,2) + pow(x[1]-0.35, 2) + pow(x[2]-0.35,2));
			return 0.15-r;
			
		}
		
		void _ic() 
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			char fname[20];
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
							
							block(ix,iy,iz).rho = _ic_func4(x);
						}
			}
		}

		
		//ADDED FROM DIEGO's REINIT
		class Integrate
		{
		public:
			double dt;
			Integrate(double dt_): dt(dt_) {}
			Integrate(const Integrate& i): dt(i.dt){}
			
			template <typename BlockType>
			inline void operator() (BlockInfo& info, BlockType& b) const
			{
				typedef BlockType B;
				
				const int n = B::sizeZ*B::sizeY*B::sizeX;
				
				//LevelsetPoint* ptrE = &(b(0));
				PLS* ptrE = &(b(0));
				for(int iE=0; iE<n; iE++, ptrE++)
				{
					//ptrE->phi += dt*ptrE->dphidt;
					ptrE->rho += dt*ptrE->drho_dt;
					ptrE->drho_dt = 0;
				}		
			}
		};
		
		bool reinit_run()
		{

			double time=0.0;
			double currTime;
			
			//ComputeRHS<nDim> computeRhs(currTime, t);
			ComputeRHS_weno computeRhs;
			Integrate integrate(5e-4);
			
			//BoundaryInfo* boundaryInfo = grid.createBoundaryInfo(computeRhs.stencil_start, computeRhs.stencil_end);
			
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
			for(int i=0; i<100; i++)
			{
				BlockProcessing::process<BlockLab>(vInfo, grid.getBlockCollection(), grid.getBoundaryInfo(), computeRhs);
				BlockProcessing::process(vInfo, grid.getBlockCollection(),  integrate);
				time += integrate.dt;
			}
			
			cout<<"Done re-init"<<endl;
			
			//delete boundaryInfo;
			
						
			return true;
		}
		
		class ComputeRHS_weno
		{
			inline const Real phi_WENO(const Real a, const Real b, const Real c, const Real d) const
			{
				const Real eps = 1e-6;
				
				const Real IS[3] = {
					13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b),
					13*(b-c)*(b-c) + 3*(b+c)*(b+c),
					13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d)
				};
				
				const Real alpha[3] = {
					1./((eps + IS[0])*(eps + IS[0])),
					6./((eps + IS[1])*(eps + IS[1])),
					3./((eps + IS[2])*(eps + IS[2])),
				};
				
				const Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
				
				const Real w[2] = {
					alpha[0]/sum_alpha,
					alpha[2]/sum_alpha
				};
				
				return 1./3*w[0]*(a-2*b+c) + 1./6*(w[1]-0.5)*(b-2*c+d);
			}
			
			inline const Real H_GODUNOV(const Real spacing, const Real phi, const Real u_plus, const Real u_minus, const Real v_plus, const Real v_minus, const Real w_plus, const Real w_minus) const
			{
				const Real s = phi/sqrt(phi*phi + spacing);
				
				if (phi>=0)
				{
					const Real term1 = max( -min(u_plus, (Real)0.0), max(u_minus, (Real)0.0) );
					const Real term2 = max( -min(v_plus, (Real)0.0), max(v_minus, (Real)0.0) );
					const Real term3 = max( -min(w_plus, (Real)0.0), max(w_minus, (Real)0.0) );
					
					return s*(sqrt(term1*term1 + term2*term2 + term3*term3)-1);
				}
				else
				{
					const Real term1 = max( -min(u_minus, (Real)0.0), max(u_plus, (Real)0.0) );
					const Real term2 = max( -min(v_minus, (Real)0.0), max(v_plus, (Real)0.0) );
					const Real term3 = max( -min(w_minus, (Real)0.0), max(w_plus, (Real)0.0) );
					
					return s*(sqrt(term1*term1 + term2*term2 + term3*term3)-1);
				}
			}
			
			template <int index, typename Field> inline const Real Dplus(Field& f, const int ix, const int iy, const int iz) const
			{ 
				return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).rho - f(ix, iy, iz).rho);
			}
			
			template <int index, typename Field> inline const Real Dminus(Field& f, const int ix, const int iy, const int iz) const
			{ 
				return (f(ix, iy, iz).rho - f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).rho );
			}
			
			template <int index, typename Field> inline const Real DplusDminus(Field& f, const int ix, const int iy, const int iz) const
			{ 
				return (f(ix + (int)(index==0), iy + (int)(index==1), iz + (int)(index==2)).rho +
						f(ix - (int)(index==0), iy - (int)(index==1), iz - (int)(index==2)).rho +
						-2*f(ix, iy,iz).rho );
			}
			
			template <typename Field>
			const Real mainRHS(Field& f, const Real spacing, const int ix, const int iy, const int iz) const
			{
				const Real u_common_term = 1./(12*spacing)*
				(-Dplus<0>(f, ix-2, iy, iz) + 7*Dplus<0>(f,ix-1, iy, iz) + 7*Dplus<0>(f,ix, iy, iz) -Dplus<0>(f, ix+1, iy, iz) );
				
				const Real v_common_term = 1./(12*spacing)*
				(-Dplus<1>(f, ix, iy-2, iz) + 7*Dplus<1>(f,ix, iy-1, iz) + 7*Dplus<1>(f,ix, iy, iz) -Dplus<1>(f, ix, iy+1, iz) );
				
				const Real w_common_term = 1./(12*spacing)*
				(-Dplus<2>(f, ix, iy, iz-2) + 7*Dplus<2>(f,ix, iy, iz-1) + 7*Dplus<2>(f,ix, iy, iz) -Dplus<2>(f, ix, iy, iz+1) );
				
				const Real phi_weno_x[2] = {
					phi_WENO(DplusDminus<0>(f, ix-2, iy, iz)/spacing,	DplusDminus<0>(f, ix-1, iy, iz)/spacing, 
							 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix+1, iy, iz)/spacing),
					phi_WENO(DplusDminus<0>(f, ix+2, iy, iz)/spacing,	DplusDminus<0>(f, ix+1, iy, iz)/spacing, 
							 DplusDminus<0>(f, ix, iy, iz)/spacing,		DplusDminus<0>(f, ix-1, iy, iz)/spacing)

				};
				
				const Real phi_weno_y[2] = {
					phi_WENO(DplusDminus<1>(f, ix, iy-2, iz)/spacing,	DplusDminus<1>(f, ix, iy-1, iz)/spacing, 
							 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy+1, iz)/spacing),
					phi_WENO(DplusDminus<1>(f, ix, iy+2, iz)/spacing,	DplusDminus<1>(f, ix, iy+1, iz)/spacing, 
							 DplusDminus<1>(f, ix, iy, iz)/spacing,		DplusDminus<1>(f, ix, iy-1, iz)/spacing)
				};
				
				const Real phi_weno_z[2] = {
					phi_WENO(DplusDminus<2>(f, ix, iy, iz-2)/spacing,	DplusDminus<2>(f, ix, iy, iz-1)/spacing, 
							 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz+1)/spacing),
					phi_WENO(DplusDminus<2>(f, ix, iy, iz+2)/spacing,	DplusDminus<2>(f, ix, iy, iz+1)/spacing, 
							 DplusDminus<2>(f, ix, iy, iz)/spacing,		DplusDminus<2>(f, ix, iy, iz-1)/spacing)
				};
				
				const Real u_minus = u_common_term - phi_weno_x[0];
				const Real u_plus = u_common_term + phi_weno_x[1];
				const Real v_minus = v_common_term - phi_weno_y[0];
				const Real v_plus = v_common_term + phi_weno_y[1];
				const Real w_minus = w_common_term - phi_weno_z[0];
				const Real w_plus = w_common_term + phi_weno_z[1];
				
				return -H_GODUNOV(spacing, f(ix,iy).rho_0, u_plus, u_minus, v_plus, v_minus, w_plus, w_minus);
			}
			
		public:
			int stencil_start[3], stencil_end[3];
			
			ComputeRHS_weno()
			{
				stencil_start[0] = -3;
				stencil_start[1] = -3;
				stencil_start[2] = -3;
				
				stencil_end[0] = 4;
				stencil_end[1] = 4;
				stencil_end[2] = 4;
			}
			
			ComputeRHS_weno(const ComputeRHS_weno& c)
			{
				stencil_start[0] = -3;
				stencil_start[1] = -3;
				stencil_start[2] = -3;
				
				stencil_end[0] = 4;
				stencil_end[1] = 4;
				stencil_end[2] = 4;
			}
			
			template<typename LabType, typename BlockType>
			inline void operator()(LabType& i, BlockInfo& info, BlockType& o) const
			{
				typedef BlockType B;
				typedef typename BlockType::ElementType E;
				
				const Real h = info.h[0];
				
				for(int iz=0; iz<B::sizeZ; iz++)
					for(int iy=0; iy<B::sizeY; iy++)
						for(int ix=0; ix<B::sizeX; ix++)
							o(ix,iy,iz).drho_dt = mainRHS(i, h , ix, iy, iz);
			}
			
		};









		
		void DumpData()
		{
			//1. create suitable filenames
			//2. dump the information of the blocks
			//3. dump the information about the data
			
			//1.
			profiler.getAgent("dumping").start();
			
			char buf[300], buf2[300];
			sprintf(buf, m_sToRenderFormat.data(), m_iCurrentDumpTime, "txt");
			sprintf(buf2, m_sToRenderFormat.data(), m_iCurrentDumpTime, "grid");
			vector<BlockInfo> vInfo = grid.getBlocksInfo();

			
			//2.
			{
				FILE * file = fopen(buf, "w");
				assert(file!=NULL);
				fprintf (file, "Wavelets: %s\n", typeid(W).name() /*"Wavelets_Interp2ndOrder\n"*/); 
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
				
				const int steStart[3] ={ W::HsSupport[0], W::HsSupport[0], W::HsSupport[0]};
				const int steEnd[3] ={  W::HsSupport[1], W::HsSupport[1], W::HsSupport[1]};
				
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
		
		void _drawLevelsetIntersections()
		{
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


		double computeDeltaTime() const
		{
			const double dx = 1.0/blockSize;
			const double CFL = 4.0;
			const double dt = dx/CFL;//8e-1;

			return dt;
		}

		void _RestartFromFile(string sFilePath, int stepid, float& output_time, Grid<W,B>& output_grid, int& out_currentDumpTime) const
		{
			//1. compute the actual time t
			//2. fill the grid from the restart file
			//3. fill the dump time

			//1.
			const double dt = computeDeltaTime();
			output_time = dt*stepid*nStepsPerSerialization;
	
			//2.
			IO_Native<W, B, levelset_projector> serializer;
			serializer.Read(output_grid, sFilePath);
			output_grid.getBoundaryInfo();
			printf("READ, size =%f MB (boundary size=%f MB)\n",  output_grid.getMemorySize(), output_grid.getBoundaryInfo().getMemorySize());
			printf("#BLOCKS = %d\n", output_grid.getBlocksInfo().size());

			//3.
			out_currentDumpTime = stepid;
		}

	public:
		Levelset(string sFormat="ReadyForVP/Grid_At_Time%03d.%s"):
		grid(blocksPerDimension,blocksPerDimension, b3D?blocksPerDimension:1, maxStencil), 
		refiner(resJump), compressor(resJump), profiler(),
		lab(), blockfwt(), viewer(true,true), t(0),
		stSorter(),m_sToRenderFormat(sFormat), m_sFormat("E:\\Smoke\\SERIALIZED_%03d"), m_iCurrentDumpTime(0)
		{
			grid.setCompressor(&compressor);
			grid.setRefiner(&refiner);

			
			stSorter.connect(grid);
			
			if (bRestartFromFile)
			{
				_RestartFromFile("E:\\Smoke\\SERIALIZED_116",116, t, grid, m_iCurrentDumpTime);
			}
			else if (!bPostProcessData)
			{
				_ic();
				Science::AutomaticRefinementForLevelsets(grid, blockfwt, toleranceIC, maxLevel);
				_ic();
				Science::AutomaticCompressionForLevelsets(grid, blockfwt, toleranceIC);         
				_ic();
			}
			
			//DUMP INITIAL FOR VP
			//DumpData();
			
			//printf("done with the IC\n");
			//IO_VTK< W, B, levelset_projector > vtk;
			//vtk.Write(grid, "babak0");
			//printf("done with the serialization to VTK\n");
			//vtk.Write(grid, "babak");
			//exit(0);
		}
		
		void Step()
		{
			if (!bPostProcessData)
			{
				const int nStep = 1;
				const float dt = computeDeltaTime();
				
				//soft refinement for memory reasons
				{
					int iMaxLoops = -1;
					bool bDoRefine = true;
					float refinement_factor = 0.75;
					
					const float memSize = grid.getMemorySize();
					if (memSize > 200.0)
					{
						iMaxLoops = 1;
						if (memSize>300) bDoRefine = false;
						else
						{
							const double lambda = (memSize-200)/(300 - 200);

							refinement_factor = refinement_factor*(1-lambda) + 1*lambda;
						}
					}

					if (bDoRefine)
						Science::AutomaticRefinementForLevelsets(grid, blockfwt, toleranceIC*refinement_factor, maxLevel, true, iMaxLoops);
				}

				for(int i=0; i<nStep; i++, t+=dt)
					_step(t, dt);
			

				Science::AutomaticCompressionForLevelsets(grid, blockfwt, toleranceIC);  

				//serialize the data
				static int serialize_counter = 0;
				if (serialize_counter++ % nStepsPerSerialization == 0)
				{
					IO_Native<W, B, levelset_projector> serializer;
					char buf[300];
					sprintf(buf, m_sFormat.data(), m_iCurrentDumpTime);
					serializer.Write(grid, buf);	
					m_iCurrentDumpTime++;
				}

				/*static int cc = 0;
				if (cc++>2)
					exit(0);*/
			}
			else
			{
				IO_Native<W, B, levelset_projector> serializer;

				const int nGrids = 181;
				for(int iGrid=0; iGrid<=nGrids; iGrid++)
				{
					printf("POSTPROCESS PHASE \n");
					char buf[300];
					sprintf(buf, m_sFormat.data(), iGrid);
					serializer.Read(grid, buf);
					

					//RefineOnceEverywhere(0);
					DumpData();
				}
				exit(0);
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

