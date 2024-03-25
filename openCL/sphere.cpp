/*
 *  sphere.cpp
 *  Synth2
 *
 *  Created by Abel on 2010-12-08.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "sphere.h"
#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>

#define nmp(x) ((x&0x01==1)?-1:1)
#define coindp(l,m) coeff[(l*(l+1) + 2*m)]
#define coindn(l,m) coeff[(l*(l+1) + 2*m + 1)]
/*
// Read the complex basis coefficients but store as real basis coefficients
void SphericalH::readCoeff(char *IFName) {
	int l, m, n;
	double re, im, nf;
	FILE *InFile = fopen(IFName, "r");
	if (!InFile)
		return;

	fscanf(InFile, "%d", &lmax);
	n = (lmax+1)*(lmax+2);
	coeff = new double[n];
	nf = 0.0;
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			fscanf(InFile, "%lf %lf", &re, &im);
			if (m == 0) {
				coindp(l,m) = re*sqrt(4.0*M_PI);
				coindn(l,m) = 0.0;
				nf += re;
			}
			else {
				nf += 2.0*re;
				coindp(l,m) = 2.0*sqrt(4.0*M_PI)*re;
				coindn(l,m) = 2.0*sqrt(4.0*M_PI)*((double)nmp(m))*im;
			}
		}
	}
	// Normalize such that sum of all real components = 1.0
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			coindp(l,m) /= nf;
			coindn(l,m) /= nf;
		}
	}	
}

// Spherical harmonics using stored coefficients
double SphericalH::SphereH(double th, double ph) {
	int l, m;
	double res = 0.0;
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			if (m == 0) {
				res += coindp(l,m)*gsl_sf_legendre_sphPlm(l, 0, cos(th));
			}
			else {
				res += coindp(l,m)*gsl_sf_legendre_sphPlm(l, m, cos(th))*cos((double)m*ph);
				res -= coindn(l,m)*gsl_sf_legendre_sphPlm(l, m, cos(th))*sin((double)m*ph);
			}
		}
	}
	return res;
}
*/

#define STEP_C (M_PI/180.0)
// Read and store the real basis coefficients
void SphericalH::readCoeff(const char *IFName) {
	int l, m, n;
	double re, im, nf;
	double spcur, ph, th;
	double spmax, phmax, thmax;
	double spmin, phmin, thmin;
	FILE *InFile = fopen(IFName, "r");
	if (!InFile)
		return;
	
	fscanf(InFile, "%d", &lmax);
	n = (lmax+1)*(lmax+2);
	coeff = new double[n];
	nf = 0.0;
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			fscanf(InFile, "%lf %lf", &re, &im);
			if (m == 0) {
				coindp(l,m) = re*sqrt(4.0*M_PI);
				coindn(l,m) = 0.0;
				nf += re;
			}
			else {
				nf += sqrt(2.0)*(re+im);
				coindp(l,m) = 2.0*sqrt(4.0*M_PI)*re;
				coindn(l,m) = 2.0*sqrt(4.0*M_PI)*im;
			}
		}
	}
	
	spmax = -1e100;
	spmin = 1e100;
	phmin = 0.0;
	phmax = 0.0;
	thmin = 0.0;
	thmax = 0.0;
	for (ph = -M_PI; ph <= M_PI; ph += STEP_C) {
		for (th = -M_PI/2.0; th <= M_PI/2.0; th += STEP_C) {
			spcur = SphereH(th, ph);
			if (spmax < spcur) {
				spmax = spcur;
				thmax = th;
				phmax = ph;
			}
			if (spmin > spcur) {
				spmin = spcur;
				thmin = th;
				phmin = ph;
			}
		}
	}
	printf("max = (%f %f %f)\nmin = (%f %f %f)\n", spmax, thmax, phmax, spmin, thmin, phmin);
	
	// Normalize such that sum of all components = 1.0
	for (l = 0; l <= lmax; l++) {
		for (m = 0; m <= l; m++) {
			coindp(l,m) /= nf;
			coindn(l,m) /= nf;
		}
	}	
}

// Spherical harmonics using stored coefficients
double SphericalH::SphereH(double th, double ph) {
	int l, m;
	double res = 0.0;
	double *spc;
	spc = new double[lmax+1];
	gsl_sf_legendre_sphPlm_array(lmax, 0, cos(th), spc);
	for (l = 0; l <= lmax; l++) {
		res += coindp(l,0)*spc[l];
	}
	for (m = 1; m <= lmax; m++) {
		gsl_sf_legendre_sphPlm_array(lmax, m, cos(th), spc);
		for (l = m; l <= lmax; l++) {
			res += spc[l-m]*(coindp(l,m)*cos((double)m*ph) + coindn(l,m)*sin((double)m*ph));
		}
	}
	delete [] spc;
	return res;
}

SphericalH::SphericalH() {
	coeff = 0x0;
}

SphericalH::~SphericalH() {
	delete [] coeff;
}
