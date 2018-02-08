/*
 *  mainlevelset.cpp
 *  MRAG
 *
 *  Created by Diego Rossinelli on 9/5/08.
 *  Remodeled by Michael Bergdorf on 9/9/08
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "GLUT/glut.h"
#include "../GridViewer.h"
#include "Levelset.h"
#include "screenshot.h"

Levelset * levelset;

static void idle(void)
{
	levelset->Step();
	glutPostRedisplay();
}

static void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	
	levelset->Render();
	
	glutSwapBuffers();
	
	static int iFrameCounter = 0;
	char buf[300];
	sprintf(buf, "img%05d.tga", iFrameCounter++);
	printf("Writing to %s\n", buf);
	gltWriteTGA(buf);
}


int main(int argc, char ** argv)
{
	const bool bVisual = true;
	
	printf("levelset\n");
	
	if(bVisual)
	{
		glutInit(&argc, argv);
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("levelset!");
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		glOrtho(-0.2, 1.2, -0.2, 1.2, -1, 1);
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
	
	levelset = new Levelset;
	
	if(bVisual)
	{
		glutMainLoop();
	}
	else
	{
		clock_t timeStart = clock();
		while (levelset->getTime()<2.0) levelset->Step();
		clock_t timeEnd = clock();
		const double t = (timeEnd - timeStart)/(double)CLOCKS_PER_SEC;
		printf("CPU Time: %f\n", t);
	}
	
	delete levelset;
	
	return 0;
}
