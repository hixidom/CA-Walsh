//g++ -o cawalsh cawalsh.cpp glut32.lib -lopengl32 -lglu32  //use this format to compile
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <windows.h>
#include "GL/glut.h"
#include "cellauto.h"

GLint winw=750,winh=750;
char current_key;
cellauto cellblock;
int t=0,initparams[]={1,200,2,64};
bool paused=false;

void resize(int x,int y){
    if(x>0)winw=x;
    if(y>0)winh=y;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,1,1,0);
}

void dispfunc() 
{
    if(!paused){
	glClear(GL_COLOR_BUFFER_BIT);		     // Clear Screen and Depth Buffer
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	cellblock.cellrecord();
	cellblock.timewalshrecord();
	
	glViewport(0,winw/2,winw/2,winh/2); cellblock.timedisp();
	glViewport(0,0,winw/2,winh/2); cellblock.groupsumdisp(3,3,0);
	glViewport(winw/2,winw/2,winw/2,winh/2); cellblock.groupproddisp(2,2,0);
	
	glViewport(winw/2,0,winw/2,winh/2); cellblock.freqdisp();
	/*
	glViewport(0,0,winw,winh);
	glColor3f(1,0,0);
	glLineWidth(2);
	glBegin(GL_LINES);
	   glVertex2f(.5,-1);
	   glVertex2f(.5,1);
	   glVertex2f(-1,.5);
	   glVertex2f(1,.5);
	glEnd();*/
	t++;
	glutSwapBuffers();
    }
}

void init() 
{
	glClearColor(0.3, 0.0, 0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,winw,winh,0);
    cellblock.init(initparams[0],initparams[1],initparams[2],initparams[3]);
}

void specialfunc(int but,int x,int y){
    if(but==GLUT_KEY_RIGHT)cellblock.rulechange(+1);
    else if(but==GLUT_KEY_LEFT)cellblock.rulechange(-1);
}

void keyfunc(unsigned char key,int x,int y){
    if(key=='p')paused=!paused;
    else if(key=='1')cellblock.reset(1);
    else if(key=='2')cellblock.reset(2);
}

int main(int argc, char **argv) 
{
    int i;
    for(i=1;i<argc;i++)initparams[i-1]=atoi(argv[i]);
    
	glutInit(&argc, argv);                                      // GLUT initialization
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);  // Display Mode
	glutInitWindowSize(winw,winh);					// set window size
	glutCreateWindow("CA Walsh");								// create Window
	glutDisplayFunc(dispfunc);									// register Display Function
	glutIdleFunc(dispfunc);						
    //glutMouseFunc(mouse);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyfunc);
	glutSpecialFunc(specialfunc);
	//glutKeyUpFunc(keyupfunc);
	init();
	glutMainLoop();												// run GLUT mainloop
	return 0;
}
