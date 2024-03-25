#ifndef __NBODY_H__
#define __NBODY_H__

#ifdef __APPLE__
#include <OpenGL/gl.h>		// Header File For The OpenGL Library
#include <GLUT/glut.h>		// Header File For The GLut Library
#include "cl.hpp"		// Header File For The OpenCL Library
#include <OpenCL/cl_gl.h>	// Header File For The CL-GL sharing Library
#else
#include <GL/glew.h>
#include <GL/gl.h>			// Header File For The OpenGL32 Library
#include <GL/glut.h>		// Header File For The GLut Library
#include <CL/cl.hpp>			// Header File For The OpenCL Library
#include <CL/cl_gl.h>		// Header File For The CL-GL sharing Library
#endif
#include "vector.h"

#define BLOCK_SIZE 64
#define PASS_SIZE 1024

#define ARRAY_HOST 0
#define ARRAY_DEVICE 1
#define ARRAY_RENDER 2


typedef struct {
	cl_uint nx;
	cl_uint nu;
	vec1 eta;
	vec1 eps2;
	vec1 dtn;
	vec1 t;
} nbody_meta;

class nbody {
private:
	vec4 reduce_v4sum(cl::Buffer &cbi, const cl_uint nri);
	vec1 reduce_v1dmin(cl::Buffer &cbi, const cl_uint nri);
	unsigned int scan_i1(const cl_uint nri);

	unsigned int fl_ArrState;
	cl_int nx;
	cl::Buffer cb_x, cb_v, cb_a, cb_J, cb_x1, cb_v1, cb_t1, cb_t2, cb_idx, cb_tmp, cb_out, cb_meta;

	cl::Kernel nbk_arr_shiftcm, nbk_arr_shiftbm;
	cl::Kernel nbk_arr_rdv4sum, nbk_arr_rdv1min;
	cl::Kernel nbk_nbody_accel1, nbk_nbody_estep, nbk_nbody_f1pred, nbk_nbody_accel2;
	cl::Kernel nbk_nbody_list, nbk_nbody_energy;
	
	cl::Context nb_context;
	cl::Program nb_prog;
	
	cl::CommandQueue nb_queue;
	size_t gbs; // global block size
	size_t lbs; // local block size
	size_t nps; // pass size
public:
	GLuint vbo;
	vec4 *x, *v, *a, *J;
	vec4 *x1, *v1;
	nbody_meta nbm;
	vec1 tc;
	int nsc;
//	size_t n; // simulation size
	
	void ArrayToDevice();
	void ArrayToHost();
	void ArrayToRender();

	void dump();
	void dump_particles();
	void s_step();
	void init();
	void initArray(vec4 *xm, vec4 *v);
	void recenter();
	vec1 time() {return nbm.t/tc; }
	void timestep();
	void timestep(unsigned int ns);
	
	static unsigned int getn(unsigned int nt);
	void init_data();
	nbody(unsigned int nt);
	~nbody();
};

#endif //__NBODY_H__
