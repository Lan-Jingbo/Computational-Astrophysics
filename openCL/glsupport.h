/*
 *  glsupport.h
 *  nbcl
 *
 *  Created by Abel on 2011-02-02.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __GLSUPPORT_H__
#define __GLSUPPORT_H__
#ifdef __APPLE__
#include <OpenGL/gl.h>		// Header File For The OpenGL Library
#include <GLUT/glut.h>		// Header File For The GLut Library
#include <OpenCL/cl.h>		// Header File For The OpenCL Library
#include <OpenCL/cl_gl.h>	// Header File For The CL-GL sharing Library
#else
#include <GL/glew.h>
// #include <GL/gl.h>			// Header File For The OpenGL32 Library
#include <GL/glut.h>		// Header File For The GLut Library
#include <CL/cl.h>		// Header File For The OpenCL Library
#include <CL/cl_gl.h>	// Header File For The CL-GL sharing Library
#endif

// void vis_GLReshape(int x, int y);
// void vis_GLKeyDown(unsigned char key, int x, int y);
// void vis_GLSpecialDown(int key, int x, int y);
extern void vis_GLSetup(int argc, char *argv[]);
extern void vis_GLAssign();
extern void s_init_counter();

#endif // __GLSUPPORT_H__
