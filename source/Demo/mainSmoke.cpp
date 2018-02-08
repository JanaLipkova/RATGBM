/*
 *  mainSmoke.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/5/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *    fka dfka
 */
#include <stdlib.h>
#include "MRAGcore/MRAGEnvironment.h"
#if _MRAG_OS == _MRAG_OS_APPLE
#include "GLUT/glut.h"
#elif _MRAG_OS == _MRAG_OS_WINDOWS
#include "GL/glew.h"
#include "GL/glut.h"
#define _USE_MATH_DEFINES
#endif

#include "../MRAGvisual/GridViewer.h"
#include "Smoke.h"
#include "screenshot.h"
GridViewer * viewer;
Smoke * smoke;  


static void idle(void)
{
	smoke->Step();
	glutPostRedisplay();
}

static void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	
	smoke->Render();
	if(smoke->getTime()>20.0) exit(0);
	//viewer->drawSketch(smoke->grid,false);
	
	glutSwapBuffers();
	
	static int iFrameCounter = 0;
	char buf[300];
	//sprintf(buf, "img%05d.tga", iFrameCounter++);
	//printf("Writing to %s\n", buf);
	//gltWriteTGA(buf);
}


int main(int argc, char ** argv)
{
	const bool bVisual = true;

	MRAG::Environment::setup();
	printf("Smoke\n");
	
	if(bVisual)
	{
		glutInit(&argc, argv);
		glutInitWindowSize(600,600);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("Smoke!");
		
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
	
	smoke = new Smoke;
	viewer = new GridViewer;
	
	if(bVisual)
	{
		glutMainLoop();
	}
	else
		while (smoke->getTime()<20.0) smoke->Step();
	
	delete smoke;

	return 0;
}




