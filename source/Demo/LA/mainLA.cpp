#include <stdlib.h>
#include "simulationLA.h"
#include <iostream>
#include "../../MRAGcore/MRAGEnvironment.h"

#ifdef	_MRAG_GLUT_VIZ
#ifdef _WIN32
#include "GL/glut.h"
#else
#include "GLUT/glut.h"
#endif
#endif

#ifdef	_MRAG_GLUT_VIZ
const int nMaxSteps = 100000;
int iCurrStep = 0;
static void display(void)
{

	glClear(GL_COLOR_BUFFER_BIT);
	
	simulation_render(true);
	simulation_run(1);
	
	
	const bool bStop = (iCurrStep++>=nMaxSteps);
	if (bStop) exit(0);
	
	glutSwapBuffers();
}

static void idle(void)
{
	glutPostRedisplay();
}
#endif

int main(int argc, char ** argv)
{
	srand(3290);
	
#ifdef _MRAG_GLUT_VIZ
	{
		glutInit(&argc, argv);
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH| GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("MRAG Refinement Test");
		
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
#endif			
	simulation_init();
	
#ifdef _MRAG_GLUT_VIZ			
	glutMainLoop();
#else
	const int nsteps = 4;
	for(int i=0; i<nsteps; i++)
		simulation_run(1);
#endif
}

