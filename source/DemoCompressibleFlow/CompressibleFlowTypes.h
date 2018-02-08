/*
 *  CompressibleFlowTypes.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/10/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


#include "../MRAGmultithreading/MRAGWrapperCUDA.h"

#include "MRAGcore/MRAGHelperMacros.h"

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

#ifdef _CUDA_SIDE
#define _call_type __device__
#else
#define _call_type 
#endif

struct FluidElement
{
	Real phi, dphidt;
	Real tmp;//, pad;
	
	void operator += (FluidElement t)
	{
		phi += t.phi;
		dphidt += t.dphidt;
	}
	
	operator Real()  { return phi; }

	FluidElement()
	{
		phi = 0.0;
		dphidt = 0.0;
	}
	
	FluidElement(const FluidElement& f):
		phi(f.phi), dphidt(f.dphidt)
	{
	}
	
	FluidElement(Real r)
	{
		phi = r;
		dphidt = 0.0;
	}

	_call_type inline void integrate(Real dt)
	{
		phi += dphidt*dt;
		dphidt = 0;
	}

	
	_call_type inline void _update_RHS(Real dt)
	{
		dphidt = tmp;
		phi -= dt*tmp;
	}
	
	
	_call_type Real evaluate_phi(Real dt)
	{
		return phi + dphidt*dt;
	}
}; 

inline FluidElement operator*(const FluidElement& p, Real v)
{
	FluidElement t;
	t.phi = p.phi*v;
	t.dphidt =  p.dphidt*v;

	return t;
}

const int blockSize = 32;
typedef Block<FluidElement, blockSize,blockSize,1> B;

template <typename T, int i> inline Real diffusion_projector_impl(const T&t)
{
	return (Real)t.phi;
}

make_projector(diffusion_projector, diffusion_projector_impl)

template<typename Container>
struct SetIC
{
	//float dt;
	char dummy;

	SetIC():dummy(0) {} 
	//SetIC(const SetIC<B>& s):dummy(0){}

	_call_type inline Real _ic_funcA(const float x[3])
	{
		const Real r[2] = {x[0]-0.5, x[1]-0.5};
		const Real IrI = sqrt(r[0]*r[0]+r[1]*r[1]);
		
		return cos(IrI*IrI*2*M_PI*10);
	}

	_call_type  void operator()(BlockInfo& info, Container& b,  int s[3],  int e[3]) 
	{
		const float x0 = info.origin[0] ;
		const float y0 = info.origin[1] ;
		const float h= info.h[0];

		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					const float x[3] = {
						x0+ix*h,
						y0+iy*h,
						0.
					};
					
					b(ix, iy, iz).phi = _ic_funcA(x);//(ix*0.01);//_ic_funcA(x);
					b(ix, iy, iz).dphidt = 0;
					b(ix, iy, iz).tmp = 0;
				}
	}
};
CudaWrapper(B,SetIC)

template<typename Container>
struct SimpleFiltering
{
	
	char stencil_start[3], stencil_end[3];
	
	_call_type inline void operator()(BlockInfo& info, Container & b, const int s[3], const int e[3])	
	{
		const float w[3] = {0.25,0.5,0.25};
		const int nPoints = (e[2]-s[2])*(e[1]-s[1])*(e[0]-s[0]);
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					b(ix, iy, iz).dphidt = b(ix-1, iy-1, iz).phi; 
				}
	}
	
	SimpleFiltering()
	{
		stencil_start[0] = stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	SimpleFiltering(const SimpleFiltering<B>& simplefiltering)
	{
		stencil_start[0] = stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
};
CudaWrapper_Gathering(B,SimpleFiltering)

template<typename Container>
struct IntegrateDiffusion
{
	float dt;

	IntegrateDiffusion(float dt_):dt(dt_){}
	IntegrateDiffusion(const IntegrateDiffusion<B>& s):dt(s.dt){}

	_call_type inline void operator()(BlockInfo& info, Container& b,  int s[3],  int e[3]) 
	{
		const float local_dt = dt;
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					/*FluidElement e = b(ix, iy, iz);

					e.phi += e.dphidt*local_dt;
					e.dphidt = 0;
					//e.integrate(local_dt);
					*/
					b(ix, iy, iz).integrate(local_dt);// = e;
				}
	}
};
CudaWrapper(B,IntegrateDiffusion)

template<typename Container>
struct ComputeDiffusionRHS
{
	float dt;
	char stencil_start[3], stencil_end[3];
	
	ComputeDiffusionRHS(float dt_): dt(dt_)
	{
		stencil_start[0] = stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	ComputeDiffusionRHS(const ComputeDiffusionRHS<B>& d): dt(d.dt)
	{
		stencil_start[0] = stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = stencil_end[1] = 2;
		stencil_end[2] = 1;
	}

	_call_type inline void operator()(BlockInfo& info, Container & b, const int s[3], const int e[3])	
	{
		float factor = 1.0/(info.h[0]*info.h[0]);

		for(int iy=s[1]; iy<e[1]; iy++)
		for(int ix=s[0]; ix<e[0]; ix++)
			b(ix,iy,0).tmp = 
				factor*(	
					b(ix-1,iy,0).evaluate_phi(dt) + b(ix+1,iy,0).evaluate_phi(dt) +
					b(ix,iy-1,0).evaluate_phi(dt) + b(ix,iy+1,0).evaluate_phi(dt) +  - 4.0*b(ix,iy,0).evaluate_phi(dt)
				);
	}
};
CudaWrapper_Gathering(B,ComputeDiffusionRHS)

template<typename Container>
struct UpdateDiffusionRHS
{
	float dt;

	UpdateDiffusionRHS(float dt_): dt(dt_){}
	UpdateDiffusionRHS(const UpdateDiffusionRHS<B>& s):dt(s.dt){}


	_call_type inline void operator()(BlockInfo& info, Container& b,  int s[3],  int e[3]) 
	{
		const float local_dt = dt;

		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					b(ix, iy, iz)._update_RHS(local_dt);
					/*FluidElement e =  b(ix,iy,iz);
					
					e.dphidt = e.tmp;
					e.phi -= local_dt*e.tmp;

					b(ix,iy,iz) = e;*/
				}
	}
};
CudaWrapper(B,UpdateDiffusionRHS)

template<typename Container>
struct ReplaceRho
{
	float dt;

	ReplaceRho(float dt_): dt(dt_){}
	ReplaceRho(const ReplaceRho<B>& s):dt(s.dt){}

	_call_type inline void operator()(BlockInfo& info, Container& b,  int s[3],  int e[3]) 
	{
		float local_dt = dt;
		for(int iz=s[2]; iz<e[2]; iz++)
			for(int iy=s[1]; iy<e[1]; iy++)
				for(int ix=s[0]; ix<e[0]; ix++)
				{
					b(ix, iy, iz).phi = b(ix,iy,iz).dphidt;
					b(ix, iy, iz).dphidt = 0;
				}
	}
};
CudaWrapper(B,ReplaceRho)

#undef _call_type
