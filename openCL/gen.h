#ifndef __GEN_H__
#define __GEN_H__

#include <istream>
#include "vector.h"
#include "sphere.h"

extern unsigned int init_file(vec4 **xm, vec4 **v, const char *IFName);

class gen_x {
public:
	vec4 off;
	vec4 rf;
	vec1 rs;
	virtual void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) = 0;
	gen_x() { }
	gen_x(std::istream &InFile) { }
};
class gen_m {
public:
	vec1 mf;
	virtual void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) = 0;
	gen_m() { }
	gen_m(std::istream &InFile) { }
};
class gen_v {
public:
	vec4 bm;
	vec1 pf;
	virtual void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) = 0;
	gen_v() { }
	gen_v(std::istream &InFile) { }
};

class gen_x_ring : public gen_x  {
public:
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_x_ring();
	gen_x_ring(std::istream &InFile);
};

class gen_x_uni : public gen_x  {
public:
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_x_uni();
	gen_x_uni(std::istream &InFile);
};

class gen_x_pow : public gen_x {
private:
	vec1 norm(vec1 a);
public:
	vec1 alpha;
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_x_pow();
	gen_x_pow(std::istream &InFile);
};

class gen_x_sph : public gen_x {
private:
	SphericalH sph;
public:
	void init_sph(const char *IFName);
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_x_sph();
	gen_x_sph(std::istream &InFile);
};

class gen_m_ring : public gen_m  {
public:
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_m_ring();
	gen_m_ring(std::istream &InFile);
};

class gen_m_uni : public gen_m  {
public:
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_m_uni();
	gen_m_uni(std::istream &InFile);
};

class gen_m_smf : public gen_m {
private:
	vec1 smfcdf(vec1 k, vec1 p, vec1 a, vec1 Ms);
public:
	vec1 alpha, M0, Ms;
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);	
	gen_m_smf();
	gen_m_smf(std::istream &InFile);
};

class gen_v_iso : public gen_v  {
public:
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_v_iso();
	gen_v_iso(std::istream &InFile);
};

class gen_v_ring : public gen_v  {
public:
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_v_ring();
	gen_v_ring(std::istream &InFile);
};

class gen_v_rot : public gen_v  {
public:
	vec1 sf;
	vec4 rv;
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);
	gen_v_rot();
	gen_v_rot(std::istream &InFile);
};

class gen_v_ivi : public gen_v {
public:
	vec1 psi;
	void gen(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n);	
	gen_v_ivi();
	gen_v_ivi(std::istream &InFile);
};

class gen {
public:
	gen_x *xg;
	gen_m *mg;
	gen_v *vg;
	void gendata(vec4 *xm, vec4 *v, unsigned int idx, unsigned int n) {
		xg->gen(xm, v, idx, n);
		mg->gen(xm, v, idx, n);
		vg->gen(xm, v, idx, n);
	}
	gen() {
		xg = 0x0;
		mg = 0x0;
		vg = 0x0;
	}
	~gen() {
		delete xg;
		delete mg;
		delete vg;
	}
};


#endif // __GEN_H__
