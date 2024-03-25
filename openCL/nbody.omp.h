#ifndef __NBODY_H__
#define __NBODY_H__

#ifdef __APPLE__
#include <OpenGL/gl.h>		// Header File For The OpenGL Library
#include <GLUT/glut.h>		// Header File For The GLut Library
#include <OpenCL/cl.hpp>		// Header File For The OpenCL Library
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
	unsigned int nx;
	unsigned int nu;
	vec1 eta;
	vec1 eps2;
	double dtn;
	double t;
} nbody_meta;

class nbody {
private:
	inline void nbody_f1pred(double *it1, vec4 *ix, vec4 *iv, vec4 *ia, vec4 *iJ, vec4 *ox, vec4 *ov);
	inline void nbody_f2pred(unsigned int *iidx, double *it1, double *it2, vec4 *ix, vec4 *iv, vec4 *ia, vec4 *iJ, vec4 *ix1, vec4 *iv1);
	inline void nbody_accel(const unsigned int gti, vec4 *ix, vec4 *iv, vec4 *oa, vec4 *oJ);
	unsigned int nbody_list(unsigned int *idx, double *t2);
	inline double nbody_tnext(double *it2);
	vec4 nbody_energy(vec4 *ix, vec4 *iv);
public:
	GLuint vbo;
	vec4 *x, *v, *a, *J;
	vec4 *x1, *v1;
	double *t1, *t2;
	unsigned int *idx;
// 	vec1 eps2, eta;
	nbody_meta nbm;
	double t, tc;
	int nsc;
// 	size_t n; // simulation size
// 	size_t gbs; // global block size
// 	size_t nbs; // local block size
// 	size_t nrx; // local blocks
// 	size_t fn; // EDFT max frequency
// 	size_t fbs; // EDFT local block size
	
	void ArrayToDevice();
	void ArrayToHost();
	void ArrayToRender();

	void dump();
	void dump_particles();
	void s_step();
	void init();
	void initArray(vec4 *xm, vec4 *v);
	void recenter();
	double time() {return t/tc;}
	void timestep();
	void timestep(unsigned int ns);
	
	static unsigned int getn(unsigned int nt);
	void init_data();
	nbody(unsigned int nt);
	~nbody();
};

#endif //__NBODY_H__
