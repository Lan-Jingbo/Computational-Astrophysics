#define __CL_ENABLE_EXCEPTIONS
#include "clsetup.h"
#include "globals.h"
#include "vector.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// Global flags
int gfl_GLCL = 0;

// Helper function to get OpenCL error string from constant
// Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
// *********************************************************************
const char* oclErrorString(cl_int error)
{
	static const char* errorString[] = {
		"CL_SUCCESS",
		"CL_DEVICE_NOT_FOUND",
		"CL_DEVICE_NOT_AVAILABLE",
		"CL_COMPILER_NOT_AVAILABLE",
		"CL_MEM_OBJECT_ALLOCATION_FAILURE",
		"CL_OUT_OF_RESOURCES",
		"CL_OUT_OF_HOST_MEMORY",
		"CL_PROFILING_INFO_NOT_AVAILABLE",
		"CL_MEM_COPY_OVERLAP",
		"CL_IMAGE_FORMAT_MISMATCH",
		"CL_IMAGE_FORMAT_NOT_SUPPORTED",
		"CL_BUILD_PROGRAM_FAILURE",
		"CL_MAP_FAILURE",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"",
		"CL_INVALID_VALUE",
		"CL_INVALID_DEVICE_TYPE",
		"CL_INVALID_PLATFORM",
		"CL_INVALID_DEVICE",
		"CL_INVALID_CONTEXT",
		"CL_INVALID_QUEUE_PROPERTIES",
		"CL_INVALID_COMMAND_QUEUE",
		"CL_INVALID_HOST_PTR",
		"CL_INVALID_MEM_OBJECT",
		"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
		"CL_INVALID_IMAGE_SIZE",
		"CL_INVALID_SAMPLER",
		"CL_INVALID_BINARY",
		"CL_INVALID_BUILD_OPTIONS",
		"CL_INVALID_PROGRAM",
		"CL_INVALID_PROGRAM_EXECUTABLE",
		"CL_INVALID_KERNEL_NAME",
		"CL_INVALID_KERNEL_DEFINITION",
		"CL_INVALID_KERNEL",
		"CL_INVALID_ARG_INDEX",
		"CL_INVALID_ARG_VALUE",
		"CL_INVALID_ARG_SIZE",
		"CL_INVALID_KERNEL_ARGS",
		"CL_INVALID_WORK_DIMENSION",
		"CL_INVALID_WORK_GROUP_SIZE",
		"CL_INVALID_WORK_ITEM_SIZE",
		"CL_INVALID_GLOBAL_OFFSET",
		"CL_INVALID_EVENT_WAIT_LIST",
		"CL_INVALID_EVENT",
		"CL_INVALID_OPERATION",
		"CL_INVALID_GL_OBJECT",
		"CL_INVALID_BUFFER_SIZE",
		"CL_INVALID_MIP_LEVEL",
		"CL_INVALID_GLOBAL_WORK_SIZE",
	};
	
	const int errorCount = sizeof(errorString) / sizeof(errorString[0]);
	
	const int index = -error;
	
	return (index >= 0 && index < errorCount) ? errorString[index] : "Unspecified Error";
}

cl::Program clsLoadProgramSource(cl::Context i_ctx, std::vector<cl::Device> i_dev, const char *progsrc) {
	unsigned int i;
	cl::Platform pf;
	cl::Program rprg;
	string pfs;
	stringstream opt;
	cl_ulong cc, ccmin;
	ccmin = 1000000000000;

	// Read source file
	ifstream sourceFile(progsrc, ifstream::in);
	string sourceCode(istreambuf_iterator<char>(sourceFile), (istreambuf_iterator<char>()));
	try {
		cl::Program::Sources source(1, make_pair(sourceCode.c_str(), sourceCode.length()+1));

		// Make program of the source code in the context
		rprg = cl::Program(i_ctx, source);

		i_dev[0].getInfo(CL_DEVICE_PLATFORM, &pf);
		for (i = 0; i < i_dev.size(); i++) {
			i_dev[i].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, &cc);
			if (cc < ccmin) ccmin = cc;
		}
		pf.getInfo(CL_PLATFORM_EXTENSIONS, &pfs);
		opt << "-cl-mad-enable ";
		opt << "-D CACHE_SIZE=" << ccmin/4 << " ";
#ifdef PREC_DOUBLE
		opt << "-D PREC_DOUBLE ";
#endif
		if(pfs.find("cl_nv_compiler_options") != std::string::npos)
			opt << " -cl-nv-verbose ";
		rprg.build(i_dev, opt.str().c_str());
		if(gf_verbose) {
			for (i = 0; i < i_dev.size(); i++) {
				rprg.getBuildInfo(i_dev[i], CL_PROGRAM_BUILD_LOG, &pfs);
				cout << "build log:" << endl << pfs << endl;
			}
		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		for (i = 0; i < i_dev.size(); i++) {
			rprg.getBuildInfo(i_dev[i], CL_PROGRAM_BUILD_LOG, &pfs);
			cout << "build log:" << endl << pfs << endl;
		}
		throw (err);
	}
	return rprg;
}

void clsDevInfo(cl_device_type i_devtype, cl::Platform &pf) {
	if(gf_verbose) {
		cl_bool av;
		unsigned int i, j;
		vector<cl::Device> devlist;
		vector<size_t> tmpn;
		cl_device_type devtype;
		cl_uint pfn, dim;
		cl_ulong memsize;
		string pfs;
		size_t maxblk;
		
		try {
			pf.getInfo(CL_PLATFORM_NAME, &pfs);
			pf.getDevices(i_devtype, &devlist);
			cout << "Platform " << pfs << ": " << devlist.size() << " devices" << endl;
			pf.getInfo(CL_PLATFORM_EXTENSIONS, &pfs);
			cout << "Extensions: " << pfs << endl;
			
			for (j = 0; j < devlist.size(); j++) {
				cout << "Device " << j << ": ";
				devlist[j].getInfo(CL_DEVICE_TYPE, &devtype);
				if (devtype == CL_DEVICE_TYPE_CPU)
					cout << "type CPU ";
				if (devtype == CL_DEVICE_TYPE_GPU)
					cout << "type GPU ";
				if (devtype == CL_DEVICE_TYPE_ACCELERATOR)
					cout << "type Accelerator ";
				devlist[j].getInfo(CL_DEVICE_NAME, &pfs);
				cout << pfs;
				devlist[j].getInfo(CL_DEVICE_VENDOR, &pfs);
				cout << " (" << pfs << ") - ";
				devlist[j].getInfo(CL_DEVICE_AVAILABLE, &av);
				if(av) {cout << "Device available" << endl;} else {cout << "Device NOT available" << endl;} 
				devlist[j].getInfo(CL_DEVICE_VERSION, &pfs);
				cout << "Device version: " << pfs << endl;
				devlist[j].getInfo(CL_DEVICE_EXTENSIONS, &pfs);
				cout << "Extensions: " << pfs << endl;
				devlist[j].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &pfn);
				cout << "Units: " << pfn << "\t";
				devlist[j].getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &pfn);
				cout << "maxclk: " << pfn << "\t";
				devlist[j].getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &dim);
				devlist[j].getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &tmpn);
				cout << "maxdim: " << dim << " ( ";
				for (i = 0; i < tmpn.size(); i++) {
					cout << tmpn[i] << " ";
				}
				devlist[j].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &maxblk);
				cout << ")\tmaxblk: " << maxblk << endl;
				cout << "Memory:" << endl;
				devlist[j].getInfo(CL_DEVICE_LOCAL_MEM_SIZE, &memsize);
				cout << "Local: " << (long)memsize << "b " << (long)memsize/1024 << "k " << (long)memsize/1048576 << "M" << endl;
				devlist[j].getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &memsize);
				cout << "Alloc: " << (long)memsize << "b " << (long)memsize/1024 << "k " << (long)memsize/1048576 << "M" << endl;
				devlist[j].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &memsize);
				cout << "Global: " << (long)memsize << "b " << (long)memsize/1024 << "k " << (long)memsize/1048576 << "M" << endl;
				devlist[j].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, &memsize);
				cout << "Cache: " << (long)memsize << "b " << (long)memsize/1024 << "k " << (long)memsize/1048576 << "M";
				devlist[j].getInfo(CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, &memsize);
				cout << " (Line: " << (long)memsize << "b)" << endl;
				cout << endl;
			}
		}
		catch (cl::Error err) {
			std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		}
	}
}

int clsGetDevScore(cl::Device &dev) {
	cl_uint cu, cl;
	cl_bool av = 0;
	
	dev.getInfo(CL_DEVICE_AVAILABLE, &av);
	if(!av)
		return 0; // This device is not available. Don't bother.
	dev.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &cu);
	dev.getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &cl);
	return cu*cl;
}
#ifdef __APPLE__
int clsGetGLCL(cl::Device &dev) {
	std::string extn;
	dev.getInfo(CL_DEVICE_EXTENSIONS, &extn);
	if(extn.find("cl_APPLE_gl_sharing") != std::string::npos)
		return 1;
	return 0;
}

int clsSetupGL(cl::Context &o_ctx, cl::Device &o_dev) {
	cl_context_properties cps[3];
	int hasCLGL = 0;
	std::vector<cl::Platform> pfl;
	std::vector<cl::Device> bestdev, devt;
	cl_uint i, j, devnt;
	CGLContextObj kCGLContext = CGLGetCurrentContext();
	CGLShareGroupObj kCGLShareGroup = CGLGetShareGroup(kCGLContext);
	
	// enumerate platforms
	try {
		cl::Platform::get(&pfl);
		// Check All devices and find best scoring device
		devnt = 0;
		for (i = 0; i < pfl.size(); i++) {
			clsDevInfo(CL_DEVICE_TYPE_GPU, pfl[i]);
			pfl[i].getDevices(CL_DEVICE_TYPE_GPU, &devt);
			for (j = 0; j < devt.size(); j++) {
				if (clsGetGLCL(devt[j])) {
					hasCLGL = 1;
					bestdev.pop_back();
					bestdev.push_back(devt[j]);
				}
			}
			devnt += devt.size();
		}
		if (devnt == 0) // no devices found
			return 0;
		if (!hasCLGL) { // no devices found with CL-GL sharing
			return 0;
		}
		
		cps[0] = CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE;
		cps[1] = (cl_context_properties)kCGLShareGroup;
		cps[2] = 0;
		o_ctx = cl::Context(bestdev, cps);
		o_dev = bestdev[0];
	}
	catch (cl::Error err) {
		std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		return 0;
	}
	gfl_GLCL = 1;
	return 1; // 1 device
}

#else
int clsGetGLCL(cl::Device &dev) {
	std::string extn;
	dev.getInfo(CL_DEVICE_EXTENSIONS, &extn);
	if(extn.find("cl_khr_gl_sharing") != std::string::npos)
		return 1;
	return 0;
}

int clsSetupGL(cl::Context &o_ctx, cl::Device &o_dev) {
	cl_context_properties cps[7];
	int hasCLGL = 0;
	std::vector<cl::Platform> pfl;
	std::vector<cl::Device> bestdev, devt;
	cl_uint i, j, devnt;
	
	// enumerate platforms
	try {
		cl::Platform::get(&pfl);
		// Check All devices and find best scoring device
		devnt = 0;
		for (i = 0; i < pfl.size(); i++) {
			clsDevInfo(CL_DEVICE_TYPE_GPU, pfl[i]);
			pfl[i].getDevices(CL_DEVICE_TYPE_GPU, &devt);
			for (j = 0; j < devt.size(); j++) {
				if (clsGetGLCL(devt[j])) {
					hasCLGL = 1;
					bestdev.clear();
					bestdev.push_back(devt[j]);
				}
			}
			devnt += devt.size();
		}
		if (devnt == 0) // no devices found
			return 0;
		if (!hasCLGL) { // no devices found with CL-GL sharing
			return 0;
		}
		
#ifdef __WIN32__
		cps[0] = CL_GL_CONTEXT_KHR;
		cps[1] = (cl_context_properties)wglGetCurrentContext();
		cps[2] = CL_WGL_HDC_KHR;
		cps[3] = (cl_context_properties)wglGetCurrentDC();
#else
		cps[0] = CL_GL_CONTEXT_KHR;
		cps[1] = (cl_context_properties)glXGetCurrentContext();
		cps[2] = CL_GLX_DISPLAY_KHR;
		cps[3] = (cl_context_properties)glXGetCurrentDisplay();
#endif
		cps[4] = 0;
	// 	err = clGetGLContextInfoKHR(cps, CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR, sizeof(cl_device_id), &bdev, NULL);
	// 	err = clGetDeviceInfo(bdev, CL_DEVICE_PLATFORM, sizeof(cl_platform_id), &bpfl, NULL);
	// 	cps[4] = CL_CONTEXT_PLATFORM;
	// 	cps[5] = (cl_context_properties)bpfl;
	// 	cps[6] = 0;
		o_ctx = cl::Context(bestdev, cps);
		o_dev = bestdev[0];
	}
	catch (cl::Error err) {
		std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		return 0;
	}

	gfl_GLCL = 1;
	return 1; // 1 device
}
#endif

int clsSetup1(cl_device_type devtype, cl::Context &o_ctx, cl::Device &o_dev) {
	std::vector<cl::Platform> pfl;
	std::vector<cl::Device> bestdev, devt;
	cl::Platform bestpfl;
	cl_uint i, j, devnt;
	int bscore, cscore;
	
	// enumerate platforms
	try {
		cl::Platform::get(&pfl);
		// Check All devices and find best scoring device
		devnt = 0;
		bscore = 0;
		for (i = 0; i < pfl.size(); i++) {
			clsDevInfo(devtype, pfl[i]);
			pfl[i].getDevices(devtype, &devt);
			for (j = 0; j < devt.size(); j++) {
				cscore = clsGetDevScore(devt[j]);
				if (cscore > bscore) {
					bscore = cscore;
					bestdev.clear();
					bestdev.push_back(devt[j]);
					bestpfl = pfl[i];
					devnt++;
				}
			}
		}
		if (devnt == 0) // no devices found
			return 0;

		cl_context_properties cps[3] = {
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)(bestpfl)(),
			0
		};
		o_ctx = cl::Context(bestdev, cps);
		o_dev = bestdev[0];
	}
	catch (cl::Error err) {
		std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		return 0;
	}

	return 1; // 1 device
}

int clsSetupAll(cl_device_type devtype, cl::Context &o_ctx, std::vector<cl::Device> &o_dev) {
	std::vector<cl::Platform> pfl;
	std::vector<cl::Device> bestdev, devt;
	cl::Platform bestpfl;
	cl_uint i, j, devnt;
	int bscore, cscore, tscore;
	
	// enumerate platforms
	try {
		cl::Platform::get(&pfl);
		
		// Check All devices and find best scoring device
		devnt = 0;
		bscore = 0;
		for (i = 0; i < pfl.size(); i++) {
			clsDevInfo(devtype, pfl[i]);
			devt.clear();
			pfl[i].getDevices(devtype, &devt);
			cscore = 0;
			for (j = 0; j < devt.size(); j++) {
				tscore = clsGetDevScore(devt[j]);
				if (tscore > 0) {
					cscore += tscore;
					devnt++;
				}
			}
			if (cscore > bscore) {
				bscore = cscore;
				bestpfl = pfl[i];
				bestdev = devt;
			}
		}
		if (devnt == 0) // no devices found
			return 0;

		cl_context_properties cps[3] = {
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)(bestpfl)(),
			0
		};
		o_ctx = cl::Context(bestdev, cps);
		o_dev = bestdev;
	}
	catch (cl::Error err) {
		std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		return 0;
	}
	return bestdev.size();
}

cl::Buffer clsGLCreateVBO(GLuint *vbo, cl::Context ctx, size_t sz) {
	cl::Buffer b;
	// create VBO
	glGenBuffers(1, vbo);
	glBindBuffer(GL_ARRAY_BUFFER, *vbo);
	
	// initialize buffer object
	glBufferData(GL_ARRAY_BUFFER, sz, NULL, GL_DYNAMIC_DRAW);
	try {
		if(gfl_GLCL) { // create OpenCL buffer from GL VBO
			b = cl::BufferGL(ctx, CL_MEM_READ_WRITE, *vbo);
		}
		else { // create standard OpenCL mem buffer
			b = cl::Buffer(ctx, CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR, sz);
		}
	}
	catch (cl::Error err) {
		std::cerr << "ERROR: " << err.what() << "(" << oclErrorString(err.err()) << ")" << std::endl;
		throw (err);
	}
	return b;
}
