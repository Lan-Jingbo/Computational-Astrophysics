/*
 *  sphere.h
 *  Synth2
 *
 *  Created by Abel on 2010-12-08.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __SPHERE_H__
#define __SPHERE_H__

#include <math.h>

class SphericalH {
public:
	int lmax;
	
	double *coeff;

	void readCoeff(const char *InFile);
	double SphereH(double th, double ph);

	SphericalH();
	~SphericalH();
};

#endif // __SPHERE_H__