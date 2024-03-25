// OpenCL setup utilities

#ifndef __CLSETUP_H__
#define __CLSETUP_H__

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include "cl.hpp"
#include <OpenCL/cl_gl.h>	// Header File For The CL-GL sharing Library
#else
#include <GL/glew.h>
// #include <GL/gl.h>
#include <GL/glx.h>
#include <CL/cl.hpp>
#include <CL/cl_gl.h>	// Header File For The CL-GL sharing Library
#endif

// #define CLSETUP_DEVINFO
#undef CLSETUP_DEVINFO

extern const char* oclErrorString(cl_int error);

extern int clsSetupGL(cl::Context &o_ctx, cl::Device &o_dev);
extern int clsSetup1(cl_device_type devtype, cl::Context &o_ctx, cl::Device &o_dev);
extern int clsSetupAll(cl_device_type devtype, cl::Context &o_ctx, std::vector<cl::Device> &o_dev);

extern cl::Program clsLoadProgramSource(cl::Context i_ctx, std::vector<cl::Device> i_dev, const char *progsrc);
extern cl::Buffer clsGLCreateVBO(GLuint *vbo, cl::Context ctx, size_t sz);


#endif //__CLSETUP_H__
