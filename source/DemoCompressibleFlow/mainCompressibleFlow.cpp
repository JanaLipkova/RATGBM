/*
 *  mainCompressibleFlow.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 10/8/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "MRAGcore/MRAGEnvironment.h"
#include "MRAGcore/MRAGCache.h"
#include <stdlib.h>
#include <stdio.h>

#if _MRAG_OS == _MRAG_OS_APPLE
	#include "GLUT/glut.h"
#elif _MRAG_OS == _MRAG_OS_WINDOWS
	#include "GL/glew.h"
	#include "GL/glut.h"
	#define _USE_MATH_DEFINES
#endif

#include "CompressibleFlow.h"
//#include "screenshot.h"

CompressibleFlow * flow;  

static void idle(void)
{
	flow->Step();
	glutPostRedisplay();
}

static void display(void)
{
	static int c= 0;
	if (c%4==0)
	{
	glClear(GL_COLOR_BUFFER_BIT);
	
	flow->Render();

	if(flow->getTime()>20.0) exit(0);
	
	glutSwapBuffers();
	
	static int iFrameCounter = 0;
	char buf[300];
	sprintf(buf, "img%05d.tga", iFrameCounter++);
	//printf("Writing to %s\n", buf);
	//gltWriteTGA(buf);
	}
	c++;
}


int main(int argc, char ** argv)
{
	MRAG::Environment::setup();
	
	
	const bool bVisual = true;

	printf("Compressible Flow!\n");
	
	if(bVisual)
	{
		glutInit(&argc, argv);
		glutInitWindowSize(600,600);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("Compressible Flow!");
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		glOrtho(-0.05, 1.05, -0.05, 1.05, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnable(GL_TEXTURE_2D);
		
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		
		glutDisplayFunc(display);
		glutIdleFunc(idle);
	}
	
	flow = new CompressibleFlow(argc, argv);
	
	if(bVisual)
	{
		glutMainLoop();
	}
	else
		while (flow->getTime()<20.0) flow->Step();
	
	delete flow;
	
	return 0;
}

