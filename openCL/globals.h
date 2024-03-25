/*
 *  globals.h
 *  nbcl
 *
 *  Created by Abel on 2011-02-02.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "nbody.h"
#include "vector.h"

// defined in parse_args.cpp
extern vec1 gp_eps; // epsilon
extern vec1 gp_dt; // base timestep
extern int gp_burst; // burst steps per frame
extern int gp_ntgt; // target number of particles
extern vec1 gp_term_tc; // crossing time
extern unsigned int gf_verbose; // verbose mode
extern unsigned int gp_term_ns; // max steps
extern unsigned int gp_term_mode; // termination mode - 0: manual, 1: auto, bit 2: set, bit 3: condition

extern void parse_args(int argc, char *argv[]);

// defined in main.cpp
extern nbody *gc_nbody; // main particle container
extern void init_nbody(nbody **pa);
extern void die(int code);

// defined in glsupport.cpp
extern int gvar_dumpintv;
extern int gvar_dumptype;

// defined in clsupport.cpp
extern int gfl_GLCL;

#endif // __GLOBALS_H__
