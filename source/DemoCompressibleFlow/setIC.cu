/*
 *  setIC.cu
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/9/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */


__global__ void setIC(float * a)
{
	a[blockIdx.x] = 0;
}