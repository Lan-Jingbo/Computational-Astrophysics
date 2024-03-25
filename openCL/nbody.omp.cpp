/*
 *  nbody.cpp
 *  nbcl
 *
 *  Created by Abel on 2011-02-01.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "nbody.h"
#include "clsetup.h"
#include "globals.h"
#include <string.h>
#include <iostream>
#include <iomanip>
#include <math.h>


void nbody::dump() {
	vec1 E, W, K, ps;
	vec4 tmp;
	ArrayToDevice();

	nbody_f1pred(t1, x, v, a, J, x1, v1);
	tmp = nbody_energy(x, v);
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
	dump();
}

unsigned int nbody::getn(unsigned int nt) {
	return lround(pow(2.0,floor(log2((double)nt))));
}

nbody::nbody(unsigned int nt) {
	nbm.eta = gp_dt;
	nbm.eps2 = gp_eps*gp_eps;
	nbm.nx = nt;
	t = 0.0;
	tc = 1.0/sqrt(nbm.nx);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	
	// initialize buffer object
	glBufferData(GL_ARRAY_BUFFER, nbm.nx*sizeof(vec4), NULL, GL_DYNAMIC_DRAW);

	x = new vec4[nbm.nx];
	v = new vec4[nbm.nx];
	a = new vec4[nbm.nx];
	J = new vec4[nbm.nx];
	x1 = new vec4[nbm.nx];
	v1 = new vec4[nbm.nx];
	
	t1 = new double[nbm.nx];
	t2 = new double[nbm.nx];
	idx = new unsigned int[nbm.nx];
	std::cout << "Integrator ready" << std::endl;
}

void nbody::ArrayToDevice() {
	return;
}

void nbody::ArrayToHost() {
	return;
}

void nbody::ArrayToRender() {
	glFinish();
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	void *tmpx = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	memcpy(tmpx, x1, nbm.nx*sizeof(vec4));
	glUnmapBuffer(GL_ARRAY_BUFFER); 
}

nbody::~nbody() {
	glDeleteBuffers(1, &vbo);
	delete [] x;
	delete [] v;
	delete [] a;
	delete [] J;
	delete [] x1;
	delete [] v1;
}

vec4 nbody::nbody_energy(vec4 *ix, vec4 *iv) {
	unsigned int j, i;
	vec1 dx, ek, ew, m;
	vec1 Ek, Ew;
	vec4 Eo;

	Ek = 0.0;
	Ew = 0.0;
#pragma omp parallel for private(ek, ew, m, dx, j) reduction(+:Ek, Ew)
	for (i = 0; i < nbm.nx; i++) {
		m = ix[i].s[3];
		ek = m*len23(iv[i]);
		ew = 0.0;
		for (j = 0; j < nbm.nx; j++) {
			if (i != j) {
				dx = len23(ix[i] - ix[j]);
				ew -= ix[j].s[3]/sqrt(dx + nbm.eps2);
			}
		}
		ew *= m;
		Ek += ek;
		Ew += ew;
	}
	Eo = 0.0;
	Eo.s[0] = Ek;
	Eo.s[1] = Ew;
	Eo/=2.0;
	return Eo;
}

inline void nbody::nbody_accel(const unsigned int gti, vec4 *ix, vec4 *iv, vec4 *oa, vec4 *oJ) {
	unsigned int j;
	vec4 dx, dv, xi, vi, xj, vj;
	vec1 rp2, at;

	xi = ix[gti].xyz();
	vi = iv[gti].xyz();
	*oa = 0.0;
	*oJ = 0.0;
	for (j = 0; j < nbm.nx; j++) {
		if (gti != j) {
			dx = (xi - ix[j]).xyz();
			dv = (vi - iv[j]).xyz();
			rp2 = 1.0/(len23(dx) + nbm.eps2);
			at = ix[j].s[3]*rp2*sqrt(rp2);
			*oa -= at*dx;
			*oJ -= at*(dv - dx*(vec1)(3.0*(dx%dv)*rp2));
		}
	}
}

inline void nbody::nbody_f1pred(double *t1, vec4 *ix, vec4 *iv, vec4 *ia, vec4 *iJ, vec4 *ox, vec4 *ov) {
	unsigned int gti;
	vec4 Ji, ai, vi, xi;
	vec1 m, dt;

#pragma omp parallel for private(Ji, ai, vi, xi, dt, m) shared (ix, iv, ia, iJ, t1, ox, ov)
	for (gti = 0; gti < nbm.nx; gti++) {
		Ji = iJ[gti].xyz();
		ai = ia[gti].xyz();
		vi = iv[gti].xyz();
		xi = ix[gti].xyz();
		dt = nbm.dtn - t1[gti];
		m = ix[gti].s[3];

		xi += ((Ji*dt/6.0f + ai*0.5f)*dt + vi)*dt;
		vi += (Ji*dt*0.5f + ai)*dt;
		xi.s[3] = m;
		ox[gti] = xi;
		ov[gti] = vi.xyz();
	}
}

inline void nbody::nbody_f2pred(unsigned int *iidx, double *it1, double *it2, vec4 *ix, vec4 *iv, vec4 *ia, vec4 *iJ, vec4 *ix1, vec4 *iv1) {
	vec4 dx, dv, da, ff2, dt;
	vec4 vi, ai, Ji;
	vec4 xi1, vi1, ai1, Ji1;
	vec1 m;
	unsigned int gti, gid;
	double dti, dtt;

#pragma omp parallel for private(Ji, ai, vi, Ji1, ai1, vi1, xi1, dt, dx, dv, da, dti, dtt, ff2, m, gti) shared (ix, iv, ia, iJ, ix1, iv1, it1, it2)
	for (gid = 0; gid < nbm.nu; gid++) {
		gti = iidx[gid];
		xi1 = ix1[gti].xyz();
		vi1 = iv1[gti].xyz();
		nbody_accel(gti, ix1, iv1, &ai1, &Ji1);

		dt = nbm.dtn - it1[gti];
		dti = it2[gti] - it1[gti];
		// functionally equivalent to 
		// dx = (((Ji+Ji1)*dt/12.0 + (ai-ai1))*dt/5.0 + (vi+vi1))*dt/2.0;
		// dv = ((Ji-Ji1)*dt/6.0 + (ai+ai1))*dt/2.0;
		Ji  = iJ[gti].xyz();
		dx = (Ji+Ji1)*dt/12.0f;
		dv = (Ji-Ji1)*dt/6.0f;
		ai  = ia[gti].xyz();
		dx = (dx + (ai-ai1))*dt*0.2f;
		dv = (dv + (ai+ai1))*dt*0.5f;
		vi  = iv[gti].xyz();
		dx = (dx + (vi+vi1))*dt*0.5f;

		m = ix[gti].s[3];
		ix[gti] += dx;
		ix[gti].s[3] = m;
		iv[gti] += dv;
		ia[gti] = ai1;
		iJ[gti] = Ji1;

		da = ai - ai1;
		ff2 = 2.0f*(3.0f*da/dt + (Ji + 2.0f*Ji1))/dt;
		dtt = sqrt(nbm.eta*sqrt(len23(ai1)/len23(ff2)));
		if (dtt < dti)
			dtt = 0.5*dti;
		else if (dtt > 2.0*dti)
			dtt = 2.0*dti;
		else
			dtt = dti;
		it1[gti] = nbm.dtn;
// 		t2[gti] = tt + dtt;
		it2[gti] = nbm.dtn + pow(2.0, floor(log2(dtt)));
	}
}

inline double nbody::nbody_tnext(double *it2) {
	unsigned int i;
	double tl, tc = 1e100;
	#pragma omp parallel shared(tc) private(tl)
	{
		tl = 1e100;
		#pragma omp for
		for (i = 0; i < nbm.nx; i++) {
			if (tl > it2[i])
				tl = it2[i];
		}
		#pragma omp critical
		{ if (tc > tl) tc = tl; }
	}
	return tc;
}

void nbody::recenter() {
	vec4 c, ct;
	vec1 cx, cy, cz;
	unsigned int i;

// 	cx = 0.0; cy = 0.0; cz = 0.0;
#pragma omp parallel shared(c) private(ct)
{
	#pragma omp single
	{ cx = 0.0; cy = 0.0; cz = 0.0;}
	#pragma omp for reduction(+:cx,cy,cz)
	for (i = 0; i < nbm.nx; i++) {
		ct = x[i].s[3]*x[i].xyz();
		cx += ct.s[0];
		cy += ct.s[1];
		cz += ct.s[2];
	}
	#pragma omp barrier
	#pragma omp single
	{
		c = vec4(cx, cy, cz, 0.0);
		c /= nbm.nx;
	}
	#pragma omp for
	for (i = 0; i < nbm.nx; i++) {
		x[i] -= c;
	}
	#pragma omp single
	{ cx = 0.0; cy = 0.0; cz = 0.0;}
	#pragma omp for reduction(+:cx,cy,cz)
	for (i = 0; i < nbm.nx; i++) {
		cx += v[i].s[0];
		cy += v[i].s[1];
		cz += v[i].s[2];
	}
	#pragma omp barrier
	#pragma omp single
	{
		c = vec4(cx, cy, cz, 0.0);
		c /= nbm.nx;
	}
	#pragma omp for
	for (i = 0; i < nbm.nx; i++) {
		v[i] -= c;
	}
}
}

unsigned int nbody::nbody_list(unsigned int *idx, double *t2) {
	unsigned int j = 0, i;
#pragma omp parallel for shared(j)
	for (i = 0; i < nbm.nx; i++) {
		if (t2[i] <= nbm.dtn) {
			#pragma omp critical (nb_list)
			{ idx[j++] = i; }
		}
	}
	return j;
}

void nbody::initArray(vec4 *xmi, vec4 *vi) {
	memcpy(x, xmi, nbm.nx*sizeof(vec4));
	memcpy(x1, xmi, nbm.nx*sizeof(vec4));
	memcpy(v, vi, nbm.nx*sizeof(vec4));
}

void nbody::init() {
	unsigned int i;
	
	nbm.dtn = pow(2.0, -18.0);
	t = 0.0;
	nsc = 0;
	#pragma omp parallel for
	for (i = 0; i < nbm.nx; i++) {
		nbody_accel(i, x, v, &a[i], &J[i]);
		t1[i] = 0.0;
		t2[i] = nbm.dtn;
	}
}


void nbody::timestep() {
	nsc++;
	nbody_f1pred(t1, x, v, a, J, x1, v1);
	nbm.nu = nbody_list(idx, t2);
	nbody_f2pred(idx, t1, t2, x, v, a, J, x1, v1);
	t = nbm.dtn;
	nbm.dtn = nbody_tnext(t2);
}

void nbody::timestep(unsigned int ns) {
	for (unsigned int i = 0; i < ns; i++)
		timestep();
}
