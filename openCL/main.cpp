#ifdef __APPLE__
#include <OpenGL/gl.h>		// Header File For The OpenGL Library
#include <GLUT/glut.h>		// Header File For The GLut Library
#include <OpenCL/cl.h>		// Header File For The OpenCL Library
#include <OpenCL/cl_gl.h>	// Header File For The CL-GL sharing Library
#else
#include <GL/glew.h>
#include <GL/gl.h>			// Header File For The OpenGL32 Library
#include <GL/glut.h>		// Header File For The GLut Library
#include <CL/cl.h>		// Header File For The OpenCL Library
#include <CL/cl_gl.h>	// Header File For The CL-GL sharing Library
#endif
#include <iostream>
#include <unistd.h>
#include <sys/time.h>
#include "gen.h"
#include "vector.h"
#include "glsupport.h"
#include "nbody.h"
#include "globals.h"

struct timeval t1, t2;
nbody *gc_nbody = 0x0;
char *IFName = 0x0;

void die(int code) {
	std::cerr << "exit " << code << std::endl;
	if(gc_nbody)
		delete gc_nbody;
	vec1 dt;
	gettimeofday(&t2, NULL);
	dt = (t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)*1e-6;
	std::cout << "time: " << dt << " seconds" << std::endl;
	exit(code);
}

void init_nbody(nbody **pa) {
	unsigned int nn;
	vec4 *xm, *v;
	gen *gen_c;
	if (IFName)
		nn = init_file(&xm, &v, IFName);
	else {
		nn = gp_ntgt;
		xm = new vec4[nn];
		v = new vec4[nn];
		gen_c = new gen();
		gen_c->xg = new gen_x_uni();
		gen_c->mg = new gen_m_uni();
		gen_c->vg = new gen_v_ivi();
//		gen_c->xg = new gen_x_ring();
//		gen_c->mg = new gen_m_ring();
//		gen_c->vg = new gen_v_ring();
		gen_c->gendata(xm, v, 0, nn);
		delete gen_c;
	}
	if (*pa == 0x0)
		*pa = new nbody(nn);
	gc_nbody->initArray(xm, v);
	gc_nbody->recenter();
	gc_nbody->init();
	gc_nbody->dump();
	s_init_counter();
	delete [] xm;
	delete [] v;
}

int main (int argc, char *argv[]) {
	vec4 cx;

	vis_GLSetup(argc, argv);
	parse_args(argc, argv);
	if (argc - optind >= 1) {
		IFName = argv[optind];
	}

	init_nbody(&gc_nbody);
	
	gettimeofday(&t1, NULL);
	
	//Let GLUT get the msgs and tell us the ones we need
	vis_GLAssign();
	return 0;
}
