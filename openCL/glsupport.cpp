/*
 *  glsupport.cpp
 *  nbcl
 *
 *  Created by Abel on 2011-02-02.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "glsupport.h"
#include "vector.h"
#include "globals.h"
#include <stdio.h>

#define PAR_FMSEC 20 // frame rate(ms/frame)
#define PAR_VPDIST 4.0 // viewport distance from origin
#define PAR_VPNEAR 1.0 // viewport near limit
#define PAR_ALPHA 0.8 
#define PAR_STEPEL (M_PI/90.0) // elevation step
#define PAR_STEPAZ (M_PI/90.0) // azimuth step

// global visualizer flags
int gfl_stereo = 0;
int gfl_anim = 0;
int gfl_axes = 1;
// global visualizer data
GLfloat gvar_az = M_PI/4.0;
GLfloat gvar_el = M_PI/4.0;
int gvar_dumpintv = 5000;
int gvar_dumptype = 0;
int gco_steps;
vec1 gvar_tstop;
int gvar_nstop;
vec1 gvar_deye = 120.0;
vec1 gvar_reye = PAR_VPDIST;

// forward declarations
static void vs_drawAxis();
static void vs_drawParticles();
static void vis_GLAnimate();

// counters init
void s_init_counter() {
	gvar_tstop = gp_term_tc;
	gvar_nstop = gp_term_ns;
	gco_steps = 0;
	if (gp_term_mode & 0x01)
		gfl_anim = 1;
}

#define DEYE (PAR_VPDIST/gvar_deye)

// displays the complete scene
static void vis_GLDisplay(void) {
	GLfloat vx, vy, vz, dx, dy;
	GLfloat tx, ty, tz;
	GLfloat ux, uy, uz;
	GLfloat ex, ey, ez;

	vis_GLAnimate();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	
	glEnable(GL_LIGHTING);
	ux = cos(gvar_az) * cos(gvar_el);
	uy = sin(gvar_az) * cos(gvar_el);
	uz = sin(gvar_el);
	vx = gvar_reye * ux;
	vy = gvar_reye * uy;
	vz = gvar_reye * uz;
	tx = (gvar_reye - 2*PAR_VPDIST) * ux;
	ty = (gvar_reye - 2*PAR_VPDIST) * uy;
	tz = (gvar_reye - 2*PAR_VPDIST) * uz;
	ex = -cos(gvar_az)*sin(gvar_el);
	ey = -sin(gvar_az)*sin(gvar_el);
	ez = cos(gvar_el);
	if (gfl_stereo == 0) {
		glLoadIdentity();	//Load a new modelview matrix -> we can apply new transformations
		gluLookAt(vx, vy, vz, tx, ty, tz, ex, ey, ez);
		
		//draw the points
		glColor3f(0.9, 0.9, 0.9);
		vs_drawParticles();
		glColor3f(0.9, 0.9, 0.0);
		vs_drawAxis();
	}
	else {
		dx = -sinf(gvar_az)*DEYE;
		dy = cosf(gvar_az)*DEYE;
		glColorMask(GL_TRUE, GL_FALSE, GL_FALSE, GL_TRUE);
		glLoadIdentity();	//Load a new modelview matrix -> we can apply new transformations
		gluLookAt(vx-dx, vy-dy, vz, tx, ty, tz, ex, ey, ez);
		glColor3f(0.9, 0.9, 0.9);
		vs_drawParticles();
		glColor3f(0.9, 0.9, 0.0);
		vs_drawAxis();
		
		glClear(GL_DEPTH_BUFFER_BIT);
		glColorMask(GL_FALSE, GL_TRUE, GL_TRUE, GL_TRUE);
		glLoadIdentity();	//Load a new modelview matrix -> we can apply new transformations
		gluLookAt(vx+dx, vy+dy, vz, tx, ty, tz, ex, ey, ez);
		glColor3f(0.9, 0.9, 0.9);
		vs_drawParticles();
		glColor3f(0.9, 0.9, 0.0);
		vs_drawAxis();
		
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	}
	glDisable(GL_LIGHTING);
	// Finish rendering
	glFlush();	
	// Swap the buffers ->make the result of rendering visible
	glutSwapBuffers();
}

static void vis_GLAnimate() {
	// Set up the next timer tick (do this first)
// 	glutTimerFunc(PAR_FMSEC, vis_GLAnimate, 0);
	vec1 ttmp;
	
	if (!gfl_anim) {
		return;
	}

	gc_nbody->timestep(gp_burst);
	gco_steps += gp_burst;
	
	if (gp_term_mode & 0x02) { // termination mode set
		if (gp_term_mode & 0x04) {
			if(gco_steps >= gvar_nstop) {
				gc_nbody->dump();
				if (gp_term_mode & 0x01)
					die(0);
				else
					gfl_anim = 0;
				gvar_nstop += gp_term_ns;
			}
		}
		else {
			ttmp = gc_nbody->time();
			if(ttmp >= gvar_tstop) {
				gc_nbody->dump();
				if (gp_term_mode & 0x01)
					die(0);
				else
					gfl_anim = 0;
				gvar_tstop = gp_term_tc*(floor(ttmp/gp_term_tc) + 1.0);
			}
		}
	}
	if (gco_steps % gvar_dumpintv == 0) {
		gc_nbody->dump();
	}

	// Force a redisplay to render the new image
// 	gc_nbody->ArrayToRender();
// 	glutPostRedisplay();
}

const GLfloat axis[][3] = {
	{0.0, 0.0, -1.2}, {0.0, 0.0, 1.2}, \
	{0.0, 0.0, 1.2}, {0.0, 0.025, 1.15}, {0.0, 0.0, 1.2}, {0.0, -0.025, 1.15}, \
	{0.0, 0.0, 1.2}, {0.025, 0.0, 1.15}, {0.0, 0.0, 1.2}, {-0.025, 0.0, 1.15}, \
	{0.0, -1.2, 0.0}, {0.0, 1.2, 0.0}, \
	{0.0, 1.2, 0.0}, {0.025, 1.15, 0.0}, {0.0, 1.2, 0.0}, {-0.025, 1.15, 0.0}, \
	{0.0, 1.2, 0.0}, {0.0, 1.15, 0.025}, {0.0, 1.2, 0.0}, {0.0, 1.15, -0.025}, \
	{-1.2, 0.0, 0.0}, {1.2, 0.0, 0.0}, \
	{1.2, 0.0, 0.0}, {1.15, 0.0, 0.025}, {1.2, 0.0, 0.0}, {1.15, 0.0, -0.025}, \
	{1.2, 0.0, 0.0}, {1.15, 0.025, 0.0}, {1.2, 0.0, 0.0}, {1.15, -0.025, 0.0}, \
};

static void vs_drawAxis() {
	if (!gfl_axes)
		return;
	glVertexPointer(3, GL_FLOAT, 0, &axis);
	glDrawArrays(GL_LINES, 0, 30);
}
static void vs_drawParticles() {
	gc_nbody->ArrayToRender();
    glBindBuffer(GL_ARRAY_BUFFER, gc_nbody->vbo);
#ifdef PREC_DOUBLE
	glVertexPointer(3, GL_DOUBLE, sizeof(vec4), 0);
#else
	glVertexPointer(3, GL_FLOAT, sizeof(vec4), 0);
#endif
	glEnableClientState(GL_VERTEX_ARRAY);
	glDrawArrays(GL_POINTS, 0, gc_nbody->nbm.nx);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// called when the window's size has changed.
static void vis_GLReshape(int x, int y) {
	if (y == 0 || x == 0) return;  //Nothing is visible then, so return
	//Set a new projection matrix
	glMatrixMode(GL_PROJECTION);  
	glLoadIdentity();
	//Angle of view:40 degrees
	gluPerspective(40.0,(GLdouble)x/(GLdouble)y, PAR_VPNEAR, 4*PAR_VPDIST + PAR_VPNEAR);
	glMatrixMode(GL_MODELVIEW);
	glViewport(0,0,x,y);  //Use the whole window for rendering
	//Adjust point size to window size
	glPointSize(GLfloat(x)/200.0);
}

// handles key presses
static void vis_GLKeyDown(unsigned char key, int x, int y) {
	switch(key)
	{
		case 'q':
		case 27:	//ESC
			gc_nbody->dump();
			die(0);
			break;
		case 'e': //reset
			gvar_az = M_PI/4.0;
			gvar_el = M_PI/4.0;
			gvar_reye = PAR_VPDIST;
			break;
		case 'a': // evolution
//			gvar_tstop = floor(gvar_tstop) + 1.0;
			gfl_anim = (gfl_anim==1)?0:1;
			break;
		case 'p': // force dump
			gc_nbody->dump();
			break;
		case 'f': // force detailed dump
			gc_nbody->dump_particles();
			break;
		case 'r': //regenerate
//			gfl_anim = 0;
//			gco_steps = 0;
			init_nbody(&gc_nbody);
			break;
		case 's': // toggle stereoscopic mode
			gfl_stereo = (gfl_stereo==1)?0:1;
			break;
//		case 'm': // stereoscopic mode parameter
//			gvar_deye -= 5.0;
//			if (gvar_deye < 5.0) { gvar_deye = 5.0; }
//			printf("deye = %f\n", gvar_deye);
//			break;
//		case 'w': // stereoscopic mode parameter
//			gvar_deye += 5.0;
//			printf("deye = %f\n", gvar_deye);
//			break;
		case 'c': // stereoscopic mode parameter
			gvar_reye -= 0.05;
			break;
		case 't': // stereoscopic mode parameter
			gvar_reye += 0.05;
			break;
		case 'l': // toggle axes
			gfl_axes = (gfl_axes==1)?0:1;
			break;
		case 'x':
			gvar_az = 0.0;
			gvar_el = 0.0;
			break;
		case 'y':
			gvar_az = M_PI/2.0;
			gvar_el = 0.0;
			break;
		case 'z':
			gvar_az = 0.0;
			gvar_el = M_PI/2.0;
			break;
			
	}
	glutPostRedisplay();
}

static void vis_GLSpecialDown(int key, int x, int y) {
	switch(key)
	{
		case  GLUT_KEY_UP :// up
			gvar_el += PAR_STEPEL;
			break;
		case  GLUT_KEY_DOWN :// down
			gvar_el -= PAR_STEPEL;
			break;
		case  GLUT_KEY_LEFT :// left
			gvar_az -= PAR_STEPAZ;
			break;
		case  GLUT_KEY_RIGHT:// right
			gvar_az += PAR_STEPAZ;
			break;
	}
	while (gvar_az > M_PI)
		gvar_az -= 2.0*M_PI;
	while (gvar_az < M_PI)
		gvar_az += 2.0*M_PI;
	if (gvar_el > M_PI/2.0)
		gvar_el = M_PI/2.0;
	if (gvar_el < -M_PI/2.0)
		gvar_el = -M_PI/2.0;
	
	glutPostRedisplay();
}

void vis_GLSetup(int argc, char *argv[]) {
	GLfloat position [] = {0.0, 0.0, 0.0, 1.0};
	GLfloat ambientLight[] = {0.5, 0.5, 0.5, 1.0};
	GLfloat diffuseLight[] = {0.0, 0.0, 0.0, 0.8};
	GLfloat specularLight[] = {0.0, 0.0, 0.0, 1.0};
	// Initialize GLUT
	glutInit(&argc, argv);
	// Use doublebuffering, RGB(A)-mode and a depth buffer
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
// 	glutInitWindowSize(1024, 768);
	glutInitWindowSize(600, 600);
	// Create a window with rendering context and everything else we need
	glutCreateWindow("Galaxies");
	// Init some state variables:
#ifndef __APPLE__
	glewInit();
#endif
	
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.1,0.1,0.1,0.0);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);	//also try GL_LINE
	
	glColorMaterial(GL_FRONT, GL_AMBIENT);
	glEnable(GL_COLOR_MATERIAL);
	
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.0);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0001);
	glEnable(GL_LIGHT0);
}

void vis_GLAssign() {
	// Assign the event-handling routines
	glutDisplayFunc(vis_GLDisplay);
	glutReshapeFunc(vis_GLReshape);
	glutKeyboardFunc(vis_GLKeyDown);
	glutSpecialFunc(vis_GLSpecialDown);
	glutIdleFunc(vis_GLDisplay);
// 	glutIdleFunc(vis_GLAnimate);

	// Start the timer
// 	glutTimerFunc(PAR_FMSEC, vis_GLAnimate, 0);
	// Start the main loop
	glutMainLoop();
}