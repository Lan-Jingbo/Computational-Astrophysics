#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "gen.h"
#include "globals.h"
#include "vector.h"
#include "sphere.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

#define PAR_RMIN 0.001

unsigned long int random_seed() {
	unsigned long int seed;
	ifstream rdev("/dev/urandom", ios::in|ios::binary);
    rdev.read((char *)(&seed), sizeof(unsigned long int));
    rdev.close();  
	return seed;
}

void gen_x_ring::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a;
	unsigned int i;

	a = vec4(-0.1, 0.0, 0.0, 0.0);
	xm[0+idx] = a*rs + off;
	a = vec4(0.1, 0.0, 0.0, 0.0);
	xm[1+idx] = a*rs + off;
	a = 0.0;	
	i = 2;
	while (i < n) {
		a.s[0] = rf.s[0]*cos((i-2.0)*2.0*M_PI/(n-2.0));
		a.s[1] = rf.s[1]*sin((i-2.0)*2.0*M_PI/(n-2.0));
		a.s[2] = 0.0;
		xm[i+idx] = a*rs + off;
		i++;
	}
}

gen_x_ring::gen_x_ring() {
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
}

gen_x_ring::gen_x_ring(istream &InFile) {
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
}

void gen_m_ring::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	unsigned int i;
	for (i = 0; i < 2; i++) {
		xm[i+idx].s[3] = mf;
	}
	for (i = 2; i < n; i++) {
		xm[i+idx].s[3] = 1e-9*mf;
	}
}

gen_m_ring::gen_m_ring() {
	mf = 1.0;
}


void gen_v_ring::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a;
	unsigned int i;
	
	a = vec4(0.0, 0.1, 0.0, 0.0);
	v[0+idx] = 2.0*a*pf/sqrt(0.2) + bm;
	a = vec4(0.0, -0.1, 0.0, 0.0);
	v[1+idx] = 2.0*a*pf/sqrt(0.2) + bm;

	i = 2;
	a = 0.0;
	while (i < n) {
		v[i+idx] = a+bm;
		i++;
	}
}

gen_v_ring::gen_v_ring() {
	bm = 0.0;
	pf = 1.0;
}

gen_v_ring::gen_v_ring(istream &InFile) {
	bm = 0.0;
	pf = 1.0;
}



void gen_x_uni::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a;
	vec1 r;
	unsigned int i;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	a = 0.0;
	i = 0;
	while (i < n) {
		a.s[0] = gsl_ran_flat(ran, -rf.s[0], rf.s[0]);
		a.s[1] = gsl_ran_flat(ran, -rf.s[1], rf.s[1]);
		a.s[2] = gsl_ran_flat(ran, -rf.s[2], rf.s[2]);
		r = (a/rf).abs();
		if (r > 1.0)
			continue;
		xm[i+idx] = a*rs + off;
		i++;
	}
	gsl_rng_free(ran);
}

gen_x_uni::gen_x_uni() {
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
}

gen_x_uni::gen_x_uni(istream &InFile) {
	vec4 val4;
	string cmd;
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
	val4.s[3] = 0.0;
	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("off") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			off = val4;
		}
		else if(cmd.compare("tri") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			rf = val4;
		}
		else if(cmd.compare("rad") == 0) {
			InFile >> rs;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}


inline vec1 gen_x_pow::norm(vec1 a) {
	if (a == 1.0) {
		return -log(PAR_RMIN);
	}
	return (pow(PAR_RMIN, 1.0-a) - 1.0)/(a-1.0);
}

void gen_x_pow::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a;
	vec1 r, nf, pp, pr;
	unsigned int i;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	cout << " rf  = ( " << rf.s[0] << " " << rf.s[1] << " " << rf.s[2] << " )" << endl;
	cout << " off = ( " << off.s[0] << " " << off.s[1] << " " << off.s[2] << " )" << endl;
	cout << " rs = " << rs << " alpha = " << alpha << endl;
	a = 0.0;
	i = 0;
	nf = norm(alpha);
	while (i < n) {
		a.s[0] = gsl_ran_flat(ran, -rf.s[0], rf.s[0]);
		a.s[1] = gsl_ran_flat(ran, -rf.s[1], rf.s[1]);
		a.s[2] = gsl_ran_flat(ran, -rf.s[2], rf.s[2]);
		r = (a/rf).abs();
		if (alpha == 0.0) {
			if (r > 1.0)
				continue;
		}
		else {
			if (r > 1.0 || r/4.0 < PAR_RMIN)
				continue;
			pp = pow(r/4.0,-alpha)/nf;
			pr = gsl_rng_uniform(ran);
			if (pr > pp)
				continue;
		}
		xm[i+idx] = a*rs + off;
		i++;
	}
	gsl_rng_free(ran);
}

gen_x_pow::gen_x_pow() {
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
	alpha = 1.8;
}

gen_x_pow::gen_x_pow(istream &InFile) {
	vec4 val4;
	string cmd;
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
	alpha = 1.8;
	val4.s[3] = 0.0;

	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("off") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			off = val4;
		}
		else if(cmd.compare("tri") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			rf = val4;
		}
		else if(cmd.compare("rad") == 0) {
			InFile >> rs;
		}
		else if(cmd.compare("alpha") == 0) {
			InFile >> alpha;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}


void gen_x_sph::init_sph(const char *IFName) {
	sph.readCoeff(IFName);	
}

void gen_x_sph::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec1 r, gaz, gel, rfs;
	vec4 a;
	unsigned int i;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	a = 0.0;
	i = 0;
	while (i < n) {
		a.s[0] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[1] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[2] = gsl_ran_flat(ran, -1.0, 1.0);
		r = (a/rf).abs();
		gaz = atan2(a.s[1], a.s[0]);
		gel = acos(a.s[2]/r);
		rfs = sph.SphereH(gel, gaz);
		if (r > rfs || r > 1.0)
			continue;
		xm[i+idx] = a*rs + off;
		i++;
	}
	gsl_rng_free(ran);
}

gen_x_sph::gen_x_sph(istream &InFile) {
	vec4 val4;
	string cmd;
	off = 0.0;
	rf = 1.0;
	rs = 1.0;
	val4.s[3] = 0.0;
	
	InFile >> cmd;
	init_sph(cmd.c_str());
	
	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("off") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			off = val4;
		}
		else if(cmd.compare("rad") == 0) {
			InFile >> rs;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}
gen_x_sph::gen_x_sph() {

}



void gen_m_uni::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	unsigned int i;
	for (i = 0; i < n; i++) {
		xm[i+idx].s[3] = mf;
	}
}

gen_m_uni::gen_m_uni() {
	mf = 1.0;
}

gen_m_uni::gen_m_uni(istream &InFile) {
	string cmd;
	mf = 1.0;
	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("mul") == 0) {
			InFile >> mf;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}


vec1 gen_m_smf::smfcdf(vec1 k, vec1 p, vec1 a, vec1 Ms) {
	return p*(gsl_sf_gamma_inc(1.0+a, 1.0/Ms) - gsl_sf_gamma_inc(1.0+a, k/Ms));
}
vec1 smf(vec1 k, vec1 p, vec1 a, vec1 Ms) {
	return (p/Ms)*pow(k/Ms, a)*exp(-k/Ms);
}

// Schechter mass function : alpha, M_star/M_min, phi(calculated from normalization)
void gen_m_smf::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec1 t1, t2, ph, mtot;
	unsigned int i;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	// Normalize the mass function
	ph = 1.0/gsl_sf_gamma_inc(1.0+alpha, M0/Ms);
	
	i = 0;
	while (i < n) {
		t1 = pow(10.0, gsl_ran_flat(ran, 0.0, 12.0));
		t2 = gsl_ran_flat(ran, 0.0, 1.0);
		if (t2 > smfcdf(t1, ph, alpha, Ms/M0))
			continue;
		xm[i+idx].s[3] = t1;
		i++;
	}
	mtot = 0.0;
	for (i = 0; i < n; i++)
		mtot += xm[i+idx].s[3];
	mtot = mf/mtot;
	for (i = 0; i < n; i++)
		xm[i+idx].s[3] *= mtot;
	
	gsl_rng_free(ran);
}

gen_m_smf::gen_m_smf() {
	mf = 1.0;
	alpha = 1.0;
	M0 = 1e10;
	Ms = 1e14;
}

gen_m_smf::gen_m_smf(istream &InFile) {
	string cmd;
	
	mf = 1.0;
	InFile >> alpha >> M0 >> Ms;
	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("mul") == 0) {
			InFile >> mf;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}



void gen_v_iso::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a;
	vec1 r, mtot;
	unsigned int i;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	i = 0;
	a = 0.0;
	while (i < n) {
		a.s[0] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[1] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[2] = gsl_ran_flat(ran, -1.0, 1.0);
		r = a.abs();
		if (r > 1.0)
			continue;
		v[i+idx] = a;
		i++;
	}
	// Remove and reinsert bulk motion
	a = 0.0;
	mtot = 0.0;
	for (i = 0; i < n; i++) {
		a += v[i+idx];
		mtot += xm[i+idx].s[3];
	}
	a /= (vec1)n;
	for (i = 0; i < n; i++) {
		v[i+idx] -= a;
		v[i+idx] *= pf*sqrt(mtot);
		v[i+idx] += bm;
	}
	gsl_rng_free(ran);
}

gen_v_iso::gen_v_iso() {
	bm = 0.0;
	pf = 1.0;
}

gen_v_iso::gen_v_iso(istream &InFile) {
	vec4 val4;
	string cmd;
	bm = 0.0;
	pf = 1.0;
	val4.s[3] = 0.0;
	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("bulk") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			bm = val4;
		}
		else if(cmd.compare("vfac") == 0) {
			InFile >> pf;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}


void gen_v_rot::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a, b;
	vec1 r, mtot;
	unsigned int i;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	i = 0;
	a = 0.0;
	while (i < n) {
		a.s[0] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[1] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[2] = gsl_ran_flat(ran, -1.0, 1.0);
		r = a.abs();
		if (r > 1.0)
			continue;
		v[i+idx] = a;
		i++;
	}
	// Remove and reinsert bulk motion
	a = 0.0;
	mtot = 0.0;
	for (i = 0; i < n; i++) {
		a += v[i+idx];
		mtot += xm[i+idx].s[3];
	}
	a /= (vec1)n;
	for (i = 0; i < n; i++) {
		v[i+idx] -= a;
		r = v[i+idx].abs();
		b = (xm[i+idx]^rv).unit()*r*(1.0-sf);
		v[i+idx] = v[i+idx]*sf + b;
		v[i+idx] *= pf*sqrt(2*mtot);
		v[i+idx] += bm;
	}
	gsl_rng_free(ran);
}

gen_v_rot::gen_v_rot() {
	bm = 0.0;
	pf = 1.0;
	sf = 0.25;
	rv = vec4(0.0, 0.0, 1.0, 0.0);
}

gen_v_rot::gen_v_rot(istream &InFile) {
	vec4 val4;
	string cmd;
	bm = 0.0;
	pf = 1.0;
	sf = 0.25;
	rv = vec4(0.0, 0.0, 1.0, 0.0);
	val4.s[3] = 0.0;

	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("bulk") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			bm = val4;
		}
		else if(cmd.compare("axis") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			rv = val4;
		}
		else if(cmd.compare("vfac") == 0) {
			InFile >> pf;
		}
		else if(cmd.compare("scatter") == 0) {
			InFile >> sf;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}


void gen_v_ivi::gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
	vec4 a;
	vec1 r, Epa, Eka, pf;
	unsigned int i, j;
	gsl_rng_env_setup();
	gsl_rng *ran = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(ran, random_seed());
	
	Epa = 0.0;
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n; j++) {
			a = xm[i+idx] - xm[j+idx];
			Epa += xm[i+idx].s[3]*xm[j+idx].s[3]/a.abs();
		}
	}
	
	i = 0;
	a = 0.0;
	Eka = 0.0;
	while (i < n) {
		a.s[0] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[1] = gsl_ran_flat(ran, -1.0, 1.0);
		a.s[2] = gsl_ran_flat(ran, -1.0, 1.0);
		r = a.abs();
		if (r > 1.0)
			continue;
		v[i+idx] = a;
		i++;
	}
	// Remove bulk motion
	a = 0.0;
	for (i = 0; i < n; i++)
		a += v[i+idx];
	a /= (vec1)n;
	for (i = 0; i < n; i++) {
		v[i+idx] -= a;
		Eka += xm[i+idx].s[3]*(v[i+idx]%v[i+idx]);
	}
	pf = sqrt(Epa/(Eka*psi));
	for (i = 0; i < n; i++) {
		v[i+idx] *= pf;
		v[i+idx] += bm;
	}
	
	gsl_rng_free(ran);
}

gen_v_ivi::gen_v_ivi() {
	bm = 0.0;
	psi = 1.0;
}

gen_v_ivi::gen_v_ivi(istream &InFile) {
	vec4 val4;
	string cmd;
	bm = 0.0;
	psi = 1.0;
	val4.s[3] = 0.0;
	while (InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("bulk") == 0) {
			InFile >> val4.s[0] >> val4.s[1] >> val4.s[2];
			bm = val4;
		}
		else if(cmd.compare("psi") == 0) {
			InFile >> psi;
		}
		else {
			InFile.seekg(-cmd.size(), ios_base::cur);
			break;
		}
	}
	if (!InFile.good() && !InFile.eof())
		throw ios_base::failure("Read Failure");
}


unsigned int init_file(vec4 **xm, vec4 **v, const char *IFName) {
	unsigned int i;
// Simulation-wide variables
	unsigned int n = -1;
	vec1 eps = 0.0001;
	vec1 mtot;
// temporary variables
	string cmd, cmd2, cmd3, cmd4;
	vec4 val4;
// current context counters
	gen *gen_c;
	unsigned int oc, nc, cc, rep;
// full cluster counters
	unsigned int ns;
	vec4 *xm1, *v1;
	vec1 eps1;
// subcluster counters
	vec4 *xm2, *v2;

	ifstream InFile(IFName, ios::in);
	
	while(InFile.good()) {
		InFile >> cmd;
		if(cmd.compare("param") == 0) {
			InFile >> n >> eps;
			cout << "Total Particles: " << n << endl;
			xm2 = new vec4[n];
			v2 = new vec4[n];
		}
		else if(cmd.compare("single") == 0) {
			cout << "Single cluster " << endl;
			xm1 = new vec4[1];
			v1 = new vec4[1];
			xm1[0] = 0.0;
			xm1[0].s[3] = n;
			v1[0] = 0.0;
			eps1 = 1.0;
			mtot = n;
			ns = 1;
			oc = 0;
			cc = 0;
		}
		else if(cmd.compare("subc") == 0) {
			InFile >> ns >> eps1;
			cout << "Subclusters: " << ns << endl;
			xm1 = new vec4[ns];
			v1 = new vec4[ns];
			gen_c = new gen();
			
			while (InFile.good()) {
				InFile >> cmd2;
				if(cmd2.compare("pos") == 0) {
					InFile >> cmd3;
					if(cmd3.compare("uni") == 0)
						gen_c->xg = new gen_x_uni(InFile);
					else if(cmd3.compare("pow") == 0)
						gen_c->xg = new gen_x_pow(InFile);
					else if(cmd3.compare("sph") == 0)
						gen_c->xg = new gen_x_sph(InFile);
					gen_c->xg->off = 0.0;
					gen_c->xg->rs = 1.0;
				}
				else if(cmd2.compare("mass") == 0) {
					InFile >> cmd3;
					if(cmd3.compare("uni") == 0)
						gen_c->mg = new gen_m_uni(InFile);
					else if(cmd3.compare("smf") == 0)
						gen_c->mg = new gen_m_smf(InFile);
					gen_c->mg->mf = (vec1)n/(vec1)ns;
				}
				else if(cmd2.compare("vel") == 0) {
					InFile >> cmd3;
					if(cmd3.compare("iso") == 0)
						gen_c->vg = new gen_v_iso(InFile);
					else if(cmd3.compare("rot") == 0)
						gen_c->vg = new gen_v_rot(InFile);
					else if(cmd3.compare("ivi") == 0)
						gen_c->vg = new gen_v_ivi(InFile);
					gen_c->vg->bm = 0.0;
				}
				else {
					InFile.seekg(-cmd2.size(), ios_base::cur);
					break;
				}
			}
			gen_c->gendata(xm1, v1, 0, ns);
			oc = 0;
			cc = 0;
			mtot = 0.0;
			for (i = 0; i < ns; i++)
				mtot += xm1[i].s[3];
			delete gen_c;
		}
		else if(cmd.compare("cluster") == 0) {
			gen_c = new gen();
			rep = 1;
			while (InFile.good()) {
				InFile >> cmd2;
				if(cmd2.compare("pos") == 0) {
					InFile >> cmd3;
					if(cmd3.compare("uni") == 0)
						gen_c->xg = new gen_x_uni(InFile);
					else if(cmd3.compare("pow") == 0)
						gen_c->xg = new gen_x_pow(InFile);
					else if(cmd3.compare("sph") == 0)
						gen_c->xg = new gen_x_sph(InFile);
				}
				else if(cmd2.compare("mass") == 0) {
					InFile >> cmd3;
					if(cmd3.compare("uni") == 0)
						gen_c->mg = new gen_m_uni(InFile);
					else if(cmd3.compare("smf") == 0)
						gen_c->mg = new gen_m_smf(InFile);
				}
				else if(cmd2.compare("vel") == 0) {
					InFile >> cmd3;
					if(cmd3.compare("iso") == 0)
						gen_c->vg = new gen_v_iso(InFile);
					else if(cmd3.compare("rot") == 0)
						gen_c->vg = new gen_v_rot(InFile);
					else if(cmd3.compare("ivi") == 0)
						gen_c->vg = new gen_v_ivi(InFile);
				}
				else if(cmd2.compare("rep") == 0) {
					InFile >> rep;
				}
				else {
					InFile.seekg(-cmd2.size(), ios_base::cur);
					break;
				}
			}
			for (i = 0; i < rep; i++) {
				if (cc >= ns)
					break;
				nc = lrint(n*xm1[cc].s[3]/mtot);
				if (oc + nc > n)
					nc = n - oc;
				cout << "Subcluster " << cc << ": " << nc << " galaxies  " << xm1[cc].s[3] << endl;
				gen_c->xg->off = xm1[cc];
				gen_c->xg->off.s[3] = 0.0;
				gen_c->xg->rf = eps1*pow(ns*xm1[cc].s[3]/mtot,1.0/3.0);
				gen_c->vg->bm = v1[cc];
				gen_c->gendata(xm2, v2, oc, nc);
				oc += nc;
				cc++;
			}
			delete gen_c;
		}

		else if(cmd.compare("END") == 0) {
			cout << "-- END OF DATA --" << endl;
			break;
		}
	}
	if(InFile.fail() && ! InFile.eof()) {
		cout << "Parse error" << endl;
		exit(-1);
	}
	InFile.close();
	
	*xm = xm2;
	*v = v2;
	gp_eps = eps;
	delete [] xm1;
	delete [] v1;
	return n;
}
