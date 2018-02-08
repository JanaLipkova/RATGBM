/*
 *  CompressibleFlow.h
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/8/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "MRAGcore/MRAGBlock.h"


#pragma once
struct FluidElement; 

namespace MRAG{ 
	template<typename W, typename B> class Grid;
	template<typename B> class BlockCollection;
	struct Wavelets_Interp4thOrder; 
	class Refiner;
}

class CompressibleFlow
{

	
	double time;
	
private:
	
public:
	CompressibleFlow(int argc, char ** argv);
	void Step();
	void Render();
	double getTime();
};
