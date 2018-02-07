
struct CPU_ComputeDiffusionRHS
{
	float dt;
	int stencil_start[3], stencil_end[3];
	
	CPU_ComputeDiffusionRHS(float dt_): dt(dt_)
	{
		stencil_start[0] = stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	CPU_ComputeDiffusionRHS(const CPU_ComputeDiffusionRHS& d): dt(d.dt)
	{
		stencil_start[0] = stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = stencil_end[1] = 2;
		stencil_end[2] = 1;
	}

	template<typename LabType, typename BlockType>
	inline void operator()(LabType& i, BlockInfo& info, BlockType& o) const
	{
		float factor = 1.0/(info.h[0]*info.h[0]);

		typedef BlockType B;
		typedef typename BlockType::ElementType E;
		
		for(int iz=0; iz<B::sizeZ; iz++)
		for(int iy=0; iy<B::sizeY; iy++)
		for(int ix=0; ix<B::sizeX; ix++)
			o(ix,iy,0).tmp = 
				factor*(	
					i(ix-1,iy,0).evaluate_phi(dt) + i(ix+1,iy,0).evaluate_phi(dt) +
					i(ix,iy-1,0).evaluate_phi(dt) + i(ix,iy+1,0).evaluate_phi(dt) +  - 4.0*i(ix,iy,0).evaluate_phi(dt)
				);
	}
};

struct CPU_IntegrateDiffusion
{
	float dt;

	CPU_IntegrateDiffusion(float dt_):dt(dt_){}

	inline void operator()(BlockInfo& info, B& b) 
	{				
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		
		FluidElement* ptrE = &(b(0));
		for(int iE=0; iE<n; iE++, ptrE++)
			ptrE->integrate(dt);
	}
};

struct CPU_UpdateDiffusionRHS
{
	float dt;

	CPU_UpdateDiffusionRHS(float dt_): dt(dt_){}

	inline void operator()(BlockInfo& info, B& b) 
	{			
		const int n = B::sizeZ*B::sizeY*B::sizeX;
		
		FluidElement* ptrE = &(b(0));
		for(int iE=0; iE<n; iE++, ptrE++)
		{
			ptrE->dphidt = ptrE->tmp;
			ptrE->phi -= dt*ptrE->tmp;
		}
	}
};


const float beta = 2.0;
struct DiffusionTestCUDA
{
	static const bool bQualityCheck = true;
	static const bool bUseGPU = true;
	static const bool bUseSpaceTimeSorter = false;
	static const bool b3D = false;
	static const int blocksPerDimension = 1;
	static const int maxLevel = 7;
	static const int resJump = 2;	

	float time;//,  dt;
	float toleranceIC;

	Grid<W, B> * grid;
	BlockCollection<B> * collection;
	Refiner refiner;
	Compressor compressor;
	SpaceTimeSorter stSorter;
	BlockProcessing_CUDA<B> processing;
	BlockFWT<W, B, diffusion_projector> blockfwt;
	Profiler profiler;

	DiffusionTestCUDA(float tol_=1e-1): 
	time(0)/*,  dt(1e-6),*/, toleranceIC(tol_),
	collection(NULL),grid(NULL)
	{
	}

	static inline Real _ic_func(Real x[3], Real time)
	{
		const Real r[2] = {x[0]-0.5, x[1]-0.5};
		const Real IrI = sqrt(r[0]*r[0]+r[1]*r[1]);
		
		return 10*cos(IrI*IrI*2*M_PI*10)*pow(sin(2*M_PI*x[0])*sin(2*M_PI*x[1]),10);
	}

	static Real _ic_test1(Real x[3],Real time)
	{
		return (exp(-8.0*M_PI*M_PI*time)*sin(2.0*M_PI*x[0])*sin(2.0*M_PI*x[1])
				+
				exp(-8.0*M_PI*M_PI*time*beta*beta)*sin(2.0*M_PI*beta*x[0])*sin(2.0*beta*M_PI*x[1]));
	}

	//template<typename Grid>
	static void _ic(MRAG::Grid<W, B>& grid)
	{
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
		
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			grid.getBlockCollection().lock(info.blockID);
			B& block = grid.getBlockCollection()[info.blockID];
			
			for(int iz=0; iz<B::sizeZ; iz++)
				for(int iy=0; iy<B::sizeY; iy++)
					for(int ix=0; ix<B::sizeX; ix++)
					{
						Real x[3];
						
						info.pos(x, ix, iy, iz);
						
						block(ix,iy,iz).phi = _ic_func(x,0);	
						//block(ix,iy,iz).phi = _ic_test1(x,0);	
					}

			grid.getBlockCollection().release(info.blockID);
		}			
	}

	void setup()
	{
		collection = NULL;

		if (bUseGPU)
			collection = new BlockCollection_CUDA<B>();


		grid = new Grid<W, B>(blocksPerDimension,blocksPerDimension,B::sizeZ>1?blocksPerDimension:1, collection);
		grid->setRefiner(&refiner);
		grid->setCompressor(&compressor);
	/*	srand(6510);
		
		*/
		IO_Native<W,B, diffusion_projector> io;

		 FILE * f = fopen("grid.txt","r");
		
		const bool bFileExists = f!=NULL;
		
		if (bFileExists && !bQualityCheck)
		{
			fclose(f);
			profiler.getAgent("Loading From File").start();
			io.Read(*grid, "grid");
			profiler.getAgent("Loading From File").stop();
		}
		else
		{
			profiler.getAgent("Creating Grid").start();
			for(int iRef = 0; iRef<(bQualityCheck?2:6); iRef++)
			{
				set<int> shouldBeRefined;
				vector<MRAG::BlockInfo> vInfo = grid->getBlocksInfo();
				
				for(int i=0; i<vInfo.size(); i++)
					if (rand()/(double)RAND_MAX > 0.4) shouldBeRefined.insert(vInfo[i].blockID);
				
				grid->refine(shouldBeRefined);
			}
			_ic(*grid);
			profiler.getAgent("Creating Grid").stop();
			/*_ic(*grid);

			Science::AutomaticRefinement< 0,0 >(*grid, blockfwt, toleranceIC, maxLevel,-1, NULL, _ic);

			_ic(*grid);
				
			Science::AutomaticCompression< 0,0 >(*grid, blockfwt, toleranceIC, -1, NULL, _ic);*/

			if (!bQualityCheck)io.Write(*grid, "grid");
		}
		stSorter.connect(*grid);
		
	}

	void _step(float time, float dt)
	{
		if (!bUseSpaceTimeSorter)
		{
			vector<BlockInfo> vInfo = grid->getBlocksInfo();

			if (bUseGPU)
			{
				SetIC<B> setIC;
		//		ComputeDiffusionRHS<B> diffusionRHS(0);
				/*UpdateDiffusionRHS<B> updateRHS(0);
				IntegrateDiffusion<B> integrate(0.000001);*/
/*static int cc = 0;
				if (cc++ == 0)processing.process_mmt(vInfo, grid->getBlockCollection(), setIC);
				processing.process_mmt(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  diffusionRHS);
				processing.process_mmt(vInfo, grid->getBlockCollection(), updateRHS);
				processing.process_mmt(vInfo, grid->getBlockCollection(), integrate);*/

				
				SimpleFiltering<B> filt;
				ReplaceRho<B> replace(0.0);
				static int cc = 0;
				if (cc++ == 0)processing.process_mmt(vInfo, grid->getBlockCollection(), setIC);
				processing.process_mmt(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  filt);
				processing.process_mmt(vInfo, grid->getBlockCollection(), replace);
			}
			else
			{
				CPU_ComputeDiffusionRHS diffusionRHS(dt);
				CPU_UpdateDiffusionRHS updateRHS(dt);
									
				BlockProcessingCPU::process<BlockLab>(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  diffusionRHS);
				BlockProcessingCPU::process(vInfo, grid->getBlockCollection(), updateRHS);

				CPU_IntegrateDiffusion integrate(dt);
				BlockProcessingCPU::process(vInfo, grid->getBlockCollection(), integrate);
			}
			//SetIC setIC;
			//SimpleFiltering<B> filt;
			//static int c = 0; 
			//if (c== 0)processing.process_mmt(vInfo, grid->getBlockCollection(), setIC);
			//c++;
			
			//processing.process_mmt(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  filt);
		/*	processing.process_mmt(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  diffusionRHS);
			processing.process_mmt(vInfo, grid->getBlockCollection(), updateRHS);
			processing.process_mmt(vInfo, grid->getBlockCollection(), integrate);*/
		}
		else
		{
			vector<BlockInfo> vInfo;
			
			double currTime, currDeltaT;
			int level;
			stSorter.startSession(dt, 4, 0);
			int c = 0;
			if (bUseGPU)
			{
				while(true)
				{
					SpaceTimeSorter::ETimeInterval type;
					const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
				
					if (type == SpaceTimeSorter::ETimeInterval_Start)
					{
						ComputeDiffusionRHS<B> diffusionRHS(currTime);
						UpdateDiffusionRHS<B> updateRHS(currTime);
						processing.process_mmt(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  diffusionRHS);
						processing.process_mmt(vInfo, grid->getBlockCollection(), updateRHS);
					}
					else
					{
						IntegrateDiffusion<B> integrate(currTime);
						processing.process_mmt(vInfo, grid->getBlockCollection(), integrate);
					}
					
					if (!bContinue) break;
				}
				
			}
			else
			{
				while(true)
				{
					SpaceTimeSorter::ETimeInterval type;
					const bool bContinue = stSorter.getBlocks(level, currDeltaT, currTime, vInfo, type);
				
					if (type == SpaceTimeSorter::ETimeInterval_Start)
					{
						CPU_ComputeDiffusionRHS diffusionRHS(currTime);
						CPU_UpdateDiffusionRHS updateRHS(currTime);
						
						BlockProcessingCPU::process<BlockLab>(vInfo, grid->getBlockCollection(), grid->getBoundaryInfo(),  diffusionRHS);
						BlockProcessingCPU::process(vInfo, grid->getBlockCollection(), updateRHS);
					}
					else
					{
						CPU_IntegrateDiffusion integrate(currTime);
						BlockProcessingCPU::process(vInfo, grid->getBlockCollection(), integrate);
					}
					
					if (!bContinue) break;
				}
			}
			stSorter.endSession();
		}
	}

	void step()
	{
		const double mu = 1;
		const int nStep = bQualityCheck?1:500;
		
		const double dx = 1.0/blockSize* (bUseSpaceTimeSorter?1:pow(2.0,-grid->getCurrentMaxLevel()));
		const double F = 1.0;//0.25;
		const float dt = F*dx*dx/(4.0*mu);//8e-1;
		
		printf("Time: %f\n",time);
		static int istep = 0;

		profiler.getAgent("Refinement Before Computing").start();
		//Science::AutomaticRefinement< 0,0 >(*grid, blockfwt, toleranceIC/2,maxLevel);
		profiler.getAgent("Refinement Before Computing").stop();

		printf("#Blocks = %d\n", grid->getBlocksInfo().size());
		_step(0, 0);
		profiler.getAgent("Computing").start();
		for(int i=0; i<nStep; i++)
		{
			printf("------------- Step %d over (%d)\n", i+1, nStep);
			_step(time, dt);
			time+=dt;
			printf("------------- END Step %d over (%d)\n", i+1, nStep);
		}
		cudaThreadSynchronize();
		profiler.getAgent("Computing").stop();

		profiler.getAgent("Compression After Computing").start();
	//	Science::AutomaticCompression< 0,0 >(*grid, blockfwt, toleranceIC);
		profiler.getAgent("Compression After Computing").stop();

		istep++;
		if(istep == 1)
		//if ( (istep>=1<<grid->getCurrentMaxLevel()+200 && bUseSpaceTimeSorter==false) || (istep == 2 && bUseSpaceTimeSorter==true))
		{
			profiler.printSummary();
			printf("End of simulation (%d)\n", 1<<grid->getCurrentMaxLevel()+1);
			printf("#Blocks = %d\n", grid->getBlocksInfo().size());
			if (!bQualityCheck) exit(0);
		//	scanf("\n");
		}
	}
};
