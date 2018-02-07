
#pragma mark UpdateRHS
struct LA_UpdateRHS_TR
{
	double dt;
	
	LA_UpdateRHS_TR(const double dt_): dt(dt_){}
	LA_UpdateRHS_TR(const LA_UpdateRHS_TR& c): dt(c.dt) {}
	
	template<typename B>
	inline void operator()(const BlockInfo& info, B& b) const
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

#pragma mark Integrate
struct LA_Integrate_TR
{
	double dt;
	
	LA_Integrate_TR(const double dt_): dt(dt_){}
	LA_Integrate_TR(const LA_Integrate_TR& c): dt(c.dt) {}
	
	template<typename B>
	inline void operator()(const BlockInfo& info, B& b) const
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


#pragma mark ComputeRHS_WENO5
struct LA_ComputeRHS_TR_WENO5
{
	int stencil_start[3], stencil_end[3];
	Real dt, sign_val;
	
	LA_ComputeRHS_TR_WENO5(const LA_ComputeRHS_TR_WENO5& c): dt(0.),sign_val(1.) 
	{
		memcpy(this, &c, sizeof(LA_ComputeRHS_TR_WENO5) );
	}
	
	LA_ComputeRHS_TR_WENO5(Real dt_, Real sign_): dt(dt_), sign_val(sign_)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
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

#pragma mark Integrate_REINIT
struct LA_Integrate_REINIT
{
	public:
		double dt;
		LA_Integrate_REINIT(double dt_): dt(dt_) {}
		LA_Integrate_REINIT(const LA_Integrate_REINIT& i): dt(i.dt){}
		
		template <typename BlockType>
		inline void operator() (const BlockInfo& info, BlockType& b) const
		{
			typedef BlockType B;
			
			const int n = B::sizeZ*B::sizeY*B::sizeX;
			
			PLS* ptrE = &(b(0));
			for(int iE=0; iE<n; iE++, ptrE++)
			{
				ptrE->rho += dt*ptrE->drho_dt;
				ptrE->drho_dt = 0;
			}		
		}
	};


#pragma mark ComputeRHS_WENO5_REINIT
struct LA_ComputeRHS_WENO5_REINIT
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
			3./((eps + IS[2])*(eps + IS[2]))
		};
		
		const Real sum_alpha = alpha[0] + alpha[1] + alpha[2];
		
		const Real w[2] = {
			alpha[0]/sum_alpha,
			alpha[2]/sum_alpha
		};
		
		return 1./3*w[0]*(a-2*b+c) + 1./6*(w[1]-0.5)*(b-2*c+d);
	}
	

	inline const double compute_sign1(const double phi, const double h, const double grad_phi_mag) const
	{
		return phi/sqrt(phi*phi + h*h);
	}
	
	inline const double compute_sign2(const double phi, const double h, const double grad_phi_mag) const
	{
		return phi/sqrt(phi*phi + h*h*grad_phi_mag*grad_phi_mag);
	}
	
	inline const double H_h(const double phi, const double h) const
	{
		if (phi<-h) return 0;
		if (phi>h) return 1;
		
		return 0.5*(1+phi/h+1/M_PI*sin(M_PI*phi/h));
	}
	
	inline const double compute_sign3(const double phi, const double h, const double grad_phi_mag) const
	{
		return 2*H_h(phi, h) - 1;
	}

	inline const Real H_GODUNOV(const Real spacing, const Real phi, const Real u_plus, const Real u_minus, const Real v_plus, const Real v_minus, const Real w_plus, const Real w_minus) const
	{
		//const Real s = phi/sqrt(phi*phi + spacing);
		
		if (phi>=0)
		{
			const Real term1 = max( -min(u_plus, (Real)0.0), max(u_minus, (Real)0.0) );
			const Real term2 = max( -min(v_plus, (Real)0.0), max(v_minus, (Real)0.0) );
			const Real term3 = max( -min(w_plus, (Real)0.0), max(w_minus, (Real)0.0) );
			const double grad_phi_mag = sqrt(term1*term1 + term2*term2 + term3*term3);
			return compute_sign1(phi, spacing, grad_phi_mag)*(sqrt(term1*term1 + term2*term2 + term3*term3)-1);
		}
		else
		{
			const Real term1 = max( -min(u_minus, (Real)0.0), max(u_plus, (Real)0.0) );
			const Real term2 = max( -min(v_minus, (Real)0.0), max(v_plus, (Real)0.0) );
			const Real term3 = max( -min(w_minus, (Real)0.0), max(w_plus, (Real)0.0) );
			const double grad_phi_mag = sqrt(term1*term1 + term2*term2 + term3*term3);
			return compute_sign1(phi, spacing, grad_phi_mag)*(sqrt(term1*term1 + term2*term2 + term3*term3)-1);
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
		
		return -H_GODUNOV(spacing, f(ix,iy,iz).rho_0, u_plus, u_minus, v_plus, v_minus, w_plus, w_minus);
	}
	
public:
	int stencil_start[3], stencil_end[3];
	
	LA_ComputeRHS_WENO5_REINIT()
	{
		stencil_start[0] = -3;
		stencil_start[1] = -3;
		stencil_start[2] = -3;
		
		stencil_end[0] = 4;
		stencil_end[1] = 4;
		stencil_end[2] = 4;
	}
	
	LA_ComputeRHS_WENO5_REINIT(const LA_ComputeRHS_WENO5_REINIT& c)
	{
		stencil_start[0] = -3;
		stencil_start[1] = -3;
		stencil_start[2] = -3;
		
		stencil_end[0] = 4;
		stencil_end[1] = 4;
		stencil_end[2] = 4;
	}
	
	template<typename LabType, typename BlockType>
	inline void operator()(LabType& i, const BlockInfo& info, BlockType& o) const
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