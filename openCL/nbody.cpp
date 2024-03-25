/*
 *  nbody.cpp
 *  nbcl
 *
 *  Created by Abel on 2011-02-01.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#define __CL_ENABLE_EXCEPTIONS
#include "nbody.h"
#include "clsetup.h"
#include "globals.h"
#include <iostream>
#include <iomanip>
#include <math.h>


void nbody::dump() {
	vec1 E, W, K, ps;
	vec4 tmp;
	ArrayToDevice();

	try {
		nbk_nbody_f1pred.setArg(0, cb_meta);
		nbk_nbody_f1pred.setArg(1, cb_t1);
		nbk_nbody_f1pred.setArg(2, cb_x);
		nbk_nbody_f1pred.setArg(3, cb_v);
		nbk_nbody_f1pred.setArg(4, cb_a);
		nbk_nbody_f1pred.setArg(5, cb_J);
		nbk_nbody_f1pred.setArg(6, cb_x1);
		nbk_nbody_f1pred.setArg(7, cb_v1);
		nb_queue.enqueueNDRangeKernel(nbk_nbody_f1pred, cl::NullRange, cl::NDRange(gbs), cl::NDRange(lbs));
		nbk_nbody_energy.setArg(0, cb_meta);
		nbk_nbody_energy.setArg(1, cb_x1);
		nbk_nbody_energy.setArg(2, cb_v1);
		nbk_nbody_energy.setArg(3, cb_tmp);
		nb_queue.enqueueNDRangeKernel(nbk_nbody_energy, cl::NullRange, cl::NDRange(gbs), cl::NDRange(lbs));
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (dump): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}
	tmp = reduce_v4sum(cb_tmp, nbm.nx);
	K = tmp.s[0];
	W = tmp.s[1];
	E = K+W;
	ps = -0.5*W/K;
	std::cout << std::showpoint << std::setprecision(5);
	std::cout << "steps = " << std::setw(8) << std::right << nsc;
	std::cout << " t = " << std::setw(8) << std::right << nbm.dtn/tc;
	std::cout << " W = " << std::setw(8) << std::right << W;
	std::cout << " K = " << std::setw(8) << std::right << K;
	std::cout << " E = " << std::setw(8) << std::right << E;
	std::cout << " psi = " << std::setw(8) << std::right << ps;
	std::cout << std::endl;
}

void nbody::dump_particles() {
// 	unsigned int i;
// 	ArrayToHost();
// 	for (i = 0; i < n; i++) {
// 		printf("%d: %6.4lg (%6.4lg %6.4lg %6.4lg)\t(%6.4lg %6.4lg %6.4lg)\t", i, x[i].w, x[i].x, x[i].y, x[i].z, v[i].x, v[i].y, v[i].z);
// 		printf(" (%6.4lg %6.4lg %6.4lg)\t%6.4lg\n", i, a[i].x, a[i].y, a[i].z, a[i].abs());
// 	}
	dump();
}

unsigned int nbody::getn(unsigned int nt) {
	return lround(pow(2.0,ceil(log2((double)nt))));
}

inline vec4 nbody::reduce_v4sum(cl::Buffer &cbi, const cl_uint nri) {
	vec4 c = 0.0;
	try {
//		if (lbs == 1) { // Probably running on CPU - do the reduction linearly
//			nb_queue.finish();
//			vec4 *rtmp = (vec4 *)nb_queue.enqueueMapBuffer(cbi, CL_TRUE, CL_MAP_READ, 0, nri*sizeof(vec4));
//			for (unsigned int j = 0; j < nri; j++) {
//				c += rtmp[j];
//			}
//			nb_queue.enqueueUnmapMemObject(cbi, rtmp);
//		}
//		else { // running on GPU
			nb_queue.finish();
			nbk_arr_rdv4sum.setArg(0, nri);
			nbk_arr_rdv4sum.setArg(1, cbi);
			nbk_arr_rdv4sum.setArg(2, cb_out);
			nbk_arr_rdv4sum.setArg(3, cl::__local(nps*sizeof(vec4)));
			nb_queue.enqueueNDRangeKernel(nbk_arr_rdv4sum, cl::NullRange, cl::NDRange(nps), cl::NDRange(nps));
			nb_queue.finish();
			nb_queue.enqueueReadBuffer(cb_out, CL_TRUE, 0, sizeof(vec4), &c);
//		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (reduce_v4sum): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}

	return c;
}

inline vec1 nbody::reduce_v1dmin(cl::Buffer &cbi, const cl_uint nri) {
	vec1 c = 1e100;
	try {
//		if (lbs == 1) { // Probably running on CPU - do the reduction linearly
//			nb_queue.finish();
//			vec1 *rtmp = (vec1 *)nb_queue.enqueueMapBuffer(cbi, CL_TRUE, CL_MAP_READ, 0, nri*sizeof(vec1));
//			for (unsigned int j = 0; j < nri; j++) {
//				if (c > rtmp[j]) {c = rtmp[j];}
//			}
//			nb_queue.enqueueUnmapMemObject(cbi, rtmp);
//		}
//		else { // running on GPU
			nb_queue.finish();
			nbk_arr_rdv1min.setArg(0, nri);
			nbk_arr_rdv1min.setArg(1, cbi);
			nbk_arr_rdv1min.setArg(2, cb_out);
			nbk_arr_rdv1min.setArg(3, cl::__local(nps*sizeof(vec1)));
			nb_queue.enqueueNDRangeKernel(nbk_arr_rdv1min, cl::NullRange, cl::NDRange(nps), cl::NDRange(nps));
			nb_queue.finish();
			nb_queue.enqueueReadBuffer(cb_out, CL_TRUE, 0, sizeof(vec1), &c);
//		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (reduce_v1dmin): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}

	return c;
}

#define arrInRender() ((fl_ArrState & 0x02) == ARRAY_RENDER)
#define arrInDevice() ((fl_ArrState & 0x01) == ARRAY_DEVICE)
#define arrInHost() ((fl_ArrState & 0x01) == ARRAY_HOST)

void nbody::ArrayToDevice() {
	std::vector<cl::Memory> mt;
	mt.push_back(cb_x1);
	if (fl_ArrState == ARRAY_DEVICE)
		return;

	try {
		if (arrInDevice() && arrInRender()) {
			// Other arrays are in the device, just get x
			nb_queue.enqueueAcquireGLObjects(&mt);
		}
		else {
			// everything is in the host
			if (gfl_GLCL) {
				ArrayToRender();
				nb_queue.enqueueAcquireGLObjects(&mt);
			}
			nb_queue.enqueueUnmapMemObject(cb_x, x);
			nb_queue.enqueueUnmapMemObject(cb_v, v);
// 			nb_queue.enqueueUnmapMemObject(cb_a, a);
// 			nb_queue.enqueueUnmapMemObject(cb_J, J);
// 			nb_queue.enqueueUnmapMemObject(cb_x1, x1);
// 			nb_queue.enqueueUnmapMemObject(cb_v1, v1);
		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (ArrayToDevice): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}

	fl_ArrState = ARRAY_DEVICE;
}

void nbody::ArrayToHost() {
	if (fl_ArrState == ARRAY_HOST)
		return;

	try {
		if (arrInHost() && arrInRender()) {
			// Other arrays are in the host, just get x1
			return;
//			glFinish();
//			glBindBuffer(GL_ARRAY_BUFFER, vbo);
//			x1 = (vec4 *)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
//			if (!x1) { std::cout << "error " << glGetError() << std::endl; die(-1);}
		}
		else {
//			if (gfl_GLCL) {
//				ArrayToRender();
//				glFinish();
//				glBindBuffer(GL_ARRAY_BUFFER, vbo);
//				x1 = (vec4 *)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
//				if (!x1) { std::cout << "error " << glGetError() << std::endl; die(-1);}
//			}
//			else {
//				x1 = (vec4 *)nb_queue.enqueueMapBuffer(cb_x, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, 0, n*sizeof(vec4));
//			}
			x = (vec4 *)nb_queue.enqueueMapBuffer(cb_x, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, 0, nbm.nx*sizeof(vec4));
			v = (vec4 *)nb_queue.enqueueMapBuffer(cb_v, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, 0, nbm.nx*sizeof(vec4));
// 			a = (vec4 *)nb_queue.enqueueMapBuffer(cb_a, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, 0, n*sizeof(vec4));
// 			J = (vec4 *)nb_queue.enqueueMapBuffer(cb_J, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, 0, n*sizeof(vec4));
// 			v1 = (vec4 *)nb_queue.enqueueMapBuffer(cb_v1, CL_TRUE, CL_MAP_READ|CL_MAP_WRITE, 0, n*sizeof(vec4));
			nb_queue.finish();
		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (ArrayToHost): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}

	fl_ArrState = ARRAY_HOST;
}

void nbody::ArrayToRender() {
	if (arrInRender())
		return;
	try {
		if (gfl_GLCL) {
			if (arrInHost()) {
				glBindBuffer(GL_ARRAY_BUFFER, vbo);
				glUnmapBuffer(GL_ARRAY_BUFFER);
			}
			if (arrInDevice()) {
				std::vector<cl::Memory> mt;
				mt.push_back(cb_x1);
				nb_queue.enqueueReleaseGLObjects(&mt);
				nb_queue.finish();
			}
			fl_ArrState = fl_ArrState | ARRAY_RENDER;
		}
		else {
			// enqueue an explicit copy, but move it into the device first
			ArrayToDevice();
			glFinish();
			glBindBuffer(GL_ARRAY_BUFFER, vbo);
			void *tmpx = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
			nb_queue.enqueueReadBuffer(cb_x1, CL_TRUE, 0, nbm.nx*sizeof(vec4), tmpx);
			nb_queue.finish();
			glUnmapBuffer(GL_ARRAY_BUFFER); 
		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (ArrayToRender): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}
}

void nbody::recenter() {
	vec4 c;

	ArrayToDevice();
	try {
		nbk_arr_shiftcm.setArg(0, cb_meta);
		nbk_arr_shiftcm.setArg(1, cb_x);
		nbk_arr_shiftcm.setArg(2, cl::__local(nps*sizeof(vec4)));
		nb_queue.enqueueNDRangeKernel(nbk_arr_shiftcm, cl::NullRange, cl::NDRange(nps), cl::NDRange(nps));

		nbk_arr_shiftbm.setArg(0, cb_meta);
		nbk_arr_shiftbm.setArg(1, cb_v);
		nbk_arr_shiftbm.setArg(2, cl::__local(nps*sizeof(vec4)));
		nb_queue.enqueueNDRangeKernel(nbk_arr_shiftbm, cl::NullRange, cl::NDRange(nps), cl::NDRange(nps));
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (recenter): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}
}

nbody::nbody(unsigned int nt) {
	std::vector<cl::Device> devl;
	cl::Device dev;
	cl_int err;
	cl_ulong maxmem, maxblk;

	nbm.nx = nt;
	nbm.eps2 = gp_eps*gp_eps;
	nbm.eta = gp_dt;
	lbs = BLOCK_SIZE;
	nps = PASS_SIZE;
	if (lbs > nbm.nx)
		lbs = getn(nbm.nx);
	gbs = lbs*(int)ceil((double)nbm.nx/(double)lbs);
	nbm.t = 0.0;
	tc = 1.0/sqrt(nbm.nx);
	nsc = 0;
	std::cout << "Initialising integrator structures:" << std::endl;
	
	// Try for GPU with sharing
	err = clsSetupGL(nb_context, dev);
	if (err == 0) { // falback to any GPU
		std::cout << "No GPU CL-GL sharing device found, falling back to any GPU" << std::endl;
		err = clsSetup1(CL_DEVICE_TYPE_GPU|CL_DEVICE_TYPE_ACCELERATOR, nb_context, dev);
	}
	if(err == 0) { // Fallback to CPU
		std::cout << "No GPU/accelerator devices found, falling back to CPU" << std::endl;
		err = clsSetup1(CL_DEVICE_TYPE_ALL, nb_context, dev);
	}
	if(err == 0) { // no devices
		std::cout << "No devices found" << std::endl;
		die(-1);
	}

	devl.push_back(dev);

	try {
		// Get device memory limits
		dev.getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &maxblk);
		dev.getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &maxmem);
		// limit particle number
		if (nbm.nx > maxmem/(40*sizeof(vec1))) {
			nbm.nx = maxmem/(40*sizeof(vec1));
			nbm.nx = lround(pow(2.0,floor(log2((double)nbm.nx))));
		}
		// limit block size
		if (lbs > maxblk) {
			lbs = lround(pow(2.0,floor(log2((double)maxblk))));
		}
		if (nps > maxblk) {
			nps = lround(pow(2.0,floor(log2((double)maxblk))));
		}
		if (nps > gbs) {
			nps = lround(pow(2.0,floor(log2((double)gbs))));
		}
		if (lbs > gbs) {
			lbs = gbs;
		}
		std::cout << "Simulation Block Size lbs = " << lbs << ", gbs = " << gbs << ", nps = " << nps << ", N = " << nbm.nx << std::endl;

		nb_queue = cl::CommandQueue(nb_context, dev, 0, &err);
		nb_prog = clsLoadProgramSource(nb_context, devl, "nbody_kern.cl");

		nbk_nbody_accel1 = cl::Kernel(nb_prog, "nbody_accel1");
		nbk_nbody_estep = cl::Kernel(nb_prog, "nbody_estep");
		nbk_nbody_f1pred = cl::Kernel(nb_prog, "nbody_f1pred");
		nbk_nbody_accel2 = cl::Kernel(nb_prog, "nbody_accel2");

		nbk_nbody_list = cl::Kernel(nb_prog, "nbody_nextlist");
		nbk_nbody_energy = cl::Kernel(nb_prog, "nbody_energy");

		nbk_arr_shiftcm = cl::Kernel(nb_prog, "array_shiftcm");
		nbk_arr_shiftbm = cl::Kernel(nb_prog, "array_shiftbm");


		nbk_arr_rdv4sum = cl::Kernel(nb_prog, "array_reduce_v4sum");
		nbk_arr_rdv1min = cl::Kernel(nb_prog, "array_reduce_v1dmin");

		std::cout << "Allocating memory" << std::endl;
		
		cb_x = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec4));
		cb_v = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec4));
		cb_a = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec4));
		cb_J = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec4));
		cb_x1 = clsGLCreateVBO(&vbo, nb_context, nbm.nx*sizeof(vec4));
		cb_v1 = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec4));

		cb_t1 = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec1));
		cb_t2 = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec1));
		cb_idx = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(cl_uint));

		cb_tmp = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, nbm.nx*sizeof(vec4));
		cb_out = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sizeof(vec4));
		cb_meta = cl::Buffer(nb_context, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sizeof(nbody_meta));
		nb_queue.enqueueWriteBuffer(cb_meta, CL_TRUE, 0, sizeof(nbody_meta), &nbm);
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (nbody): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}

	if(gfl_GLCL)
		fl_ArrState = ARRAY_RENDER|ARRAY_DEVICE;
	else
		fl_ArrState = ARRAY_DEVICE;

	std::cout << "Integrator ready" << std::endl;
}

nbody::~nbody() {
// 	ArrayToRender();
// 	clReleaseCommandQueue(nb_queue);
	glDeleteBuffers(1, &vbo);
// 	clReleaseMemObject(cb_v);
// 	clReleaseMemObject(cb_a);
// 	clReleaseMemObject(cb_m);
// 	clReleaseMemObject(cb_J);
// 	clReleaseMemObject(cb_x1);
// 	clReleaseMemObject(cb_v1);
// 	clReleaseMemObject(cb_tmp);
}

void nbody::initArray(vec4 *xm, vec4 *v) {
	ArrayToDevice();
	try {
		nb_queue.enqueueWriteBuffer(cb_x, CL_FALSE, 0, nbm.nx*sizeof(vec4), xm);
		nb_queue.enqueueWriteBuffer(cb_v, CL_FALSE, 0, nbm.nx*sizeof(vec4), v);
		nb_queue.enqueueCopyBuffer(cb_x, cb_x1, 0, 0, nbm.nx*sizeof(vec4));
		nb_queue.finish();
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (initArray): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}
}

void nbody::init() {
	ArrayToDevice();
	nsc = 0;
	nbm.t = 0.0;
	nbm.dtn = pow(2.0, -18.0);
	try {
		nb_queue.enqueueWriteBuffer(cb_meta, CL_FALSE, 0, sizeof(nbody_meta), &nbm);
		nbk_nbody_estep.setArg(0, cb_meta);
		nbk_nbody_estep.setArg(1, cb_t1);
		nbk_nbody_estep.setArg(2, cb_t2);
		nb_queue.enqueueNDRangeKernel(nbk_nbody_estep, cl::NullRange, cl::NDRange(gbs), cl::NDRange(lbs));

		nbk_nbody_accel1.setArg(0, cb_meta);
		nbk_nbody_accel1.setArg(1, cb_x);
		nbk_nbody_accel1.setArg(2, cb_v);
		nbk_nbody_accel1.setArg(3, cb_a);
		nbk_nbody_accel1.setArg(4, cb_J);
		nb_queue.enqueueNDRangeKernel(nbk_nbody_accel1, cl::NullRange, cl::NDRange(gbs), cl::NDRange(lbs));
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (init): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}
}

void nbody::timestep() {
	ArrayToDevice();
	nsc++;
	try {
		nbk_nbody_list.setArg(0, cb_meta);
		nbk_nbody_list.setArg(1, cb_t2);
		nbk_nbody_list.setArg(2, cb_idx);
		nbk_nbody_list.setArg(3, cl::__local(nps*sizeof(vec1)));
		nb_queue.enqueueNDRangeKernel(nbk_nbody_list, cl::NullRange, cl::NDRange(nps), cl::NDRange(nps));

// 		nb_queue.enqueueReadBuffer(cb_meta, CL_TRUE, 0, sizeof(nbody_meta), &nbm);
// 		cl_uint *tmp1 = (cl_uint *)nb_queue.enqueueMapBuffer(cb_idx, CL_TRUE, CL_MAP_READ, 0, nbm.nx*sizeof(cl_uint));
// 		std::cout << nbm.nu << ":";
// 		for (cl_uint j = 0; j < nbm.nx; j++)
// 			std::cout << "\t" << tmp1[j];
// 		std::cout << std::endl;
// 		nb_queue.enqueueUnmapMemObject(cb_idx, tmp1);

		nbk_nbody_f1pred.setArg(0, cb_meta);
		nbk_nbody_f1pred.setArg(1, cb_t1);
		nbk_nbody_f1pred.setArg(2, cb_x);
		nbk_nbody_f1pred.setArg(3, cb_v);
		nbk_nbody_f1pred.setArg(4, cb_a);
		nbk_nbody_f1pred.setArg(5, cb_J);
		nbk_nbody_f1pred.setArg(6, cb_x1);
		nbk_nbody_f1pred.setArg(7, cb_v1);
		nb_queue.enqueueNDRangeKernel(nbk_nbody_f1pred, cl::NullRange, cl::NDRange(gbs), cl::NDRange(lbs));

		nbk_nbody_accel2.setArg(0, cb_meta);
		nbk_nbody_accel2.setArg(1, cb_idx);
		nbk_nbody_accel2.setArg(2, cb_t1);
		nbk_nbody_accel2.setArg(3, cb_t2);
		nbk_nbody_accel2.setArg(4, cb_x);
		nbk_nbody_accel2.setArg(5, cb_v);
		nbk_nbody_accel2.setArg(6, cb_a);
		nbk_nbody_accel2.setArg(7, cb_J);
		nbk_nbody_accel2.setArg(8, cb_x1);
		nbk_nbody_accel2.setArg(9, cb_v1);
		nb_queue.enqueueNDRangeKernel(nbk_nbody_accel2, cl::NullRange, cl::NDRange(gbs), cl::NDRange(lbs));

		nb_queue.enqueueReadBuffer(cb_meta, CL_FALSE, 0, sizeof(nbody_meta), &nbm);
	}
	catch (cl::Error err) {
		std::cerr << "ERROR (timestep): " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		die(-1);
	}
}

void nbody::timestep(unsigned int ns) {
	for (unsigned int i = 0; i < ns; i++)
		timestep();
}
