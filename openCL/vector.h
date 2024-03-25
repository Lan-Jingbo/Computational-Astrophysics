// vector data structures

#ifndef __VECTOR_H__
#define __VECTOR_H__
#include <math.h>
// #define PREC_DOUBLE // Double precision
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#include <emmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif

typedef double v4df __attribute__ ((vector_size (4*sizeof(double))));
typedef float v4sf __attribute__ ((vector_size (4*sizeof(float))));

//typedef v4df vec4;
#ifdef PREC_DOUBLE
typedef double vec1;
#else
typedef float vec1;
#endif

struct vec4 {
	union {
#ifdef PREC_DOUBLE
		v4df v;
		__m128d xmm[2];
		cl_double4 c;
#else
		v4sf v;
		__m128 xmm;
		cl_float4 c;
#endif
		vec1 s[4];
//		struct { vec1 s0, s1, s2, s3; };
//		struct { vec1 x, y, z, w; };
	};
	//	inline vec4 (vec1 sf) { x = sf; y = sf; z = sf; w = sf; }
#ifdef PREC_DOUBLE
	inline vec4 () { v = (v4df){0.0, 0.0, 0.0, 0.0}; }
	inline vec4 (v4df sf) { v = sf;}
	inline vec4 (vec1 sf) { xmm[0] = _mm_set1_pd(sf); xmm[1] = _mm_set1_pd(sf); }
	inline void operator= (const vec1 *sf) { v = *(v4df*)sf;}
#else
	inline vec4 () { v = (v4sf){0.0, 0.0, 0.0, 0.0}; }
	inline vec4 (v4sf sf) { v = sf;}
	inline vec4 (vec1 sf) { v = _mm_set1_ps(sf); }
	inline void operator= (const vec1 *sf) { v = *(v4sf*)sf;}
#endif
	inline vec4 (vec1 fx, vec1 fy, vec1 fz, vec1 fw) { s[0] = fx; s[1] = fy; s[2] = fz; s[3] = fw; }
	
#ifdef __SSE3__
#ifdef PREC_DOUBLE
	inline vec1 sum () const {
		vec1 cf;
		__m128d t = _mm_add_pd(xmm[0], xmm[1]);
		t = _mm_hadd_pd(t, t);
		_mm_store_sd(&cf, t);
		return cf;
	}
	inline vec1 sum3 () const {
		vec1 cf;
		__m128d t = _mm_hadd_pd(xmm[0], xmm[0]);
		t = _mm_add_pd(t, xmm[1]);
		_mm_store_sd(&cf, t);
		return cf;
	}
#else
	inline vec1 sum () const {
		vec1 cf;
		__m128 t = _mm_hadd_ps(xmm, xmm);
		t = _mm_hadd_ps(t, t);
		_mm_store_ss(&cf, t);
		return cf;
	}
	inline vec1 sum3 () const {
		vec1 cf;
		__m128 t = _mm_hadd_ps(xmm, xmm);
		_mm_store_ss(&cf, t);
		cf += s[2];
		return cf;
	}
#endif
#else	
	inline vec1 sum () const { return s[0] + s[1] + s[2] + s[3]; }
	inline vec1 sum3 () const { return s[0] + s[1] + s[2]; }
#endif
	
	inline void operator= (const vec4 &sf) { v = sf.v; }
	
// 	inline vec4 xyz() const { vec4 a = *this; a.s[3] = 0.0; return a; }
	inline vec4 xyz() const { return vec4(s[0], s[1], s[2], 0.0); }
	inline vec1 abs() const { return sqrt((*this * *this).sum3()); }
	inline vec4 unit() const { return *this/(this->abs()); } // normalize
	inline vec1 operator% (const vec4 &sf) const { return ((*this)*sf).sum3(); } // 3D dot product
	inline vec4 operator+ (const vec4 &sf) const { return vec4(v+sf.v); }
	inline vec4 operator- (const vec4 &sf) const { return vec4(v-sf.v); }
	inline vec4 operator* (const vec4 &sf) const { return vec4(v*sf.v); }
	inline vec4 operator/ (const vec4 &sf) const { return vec4(v/sf.v); }
	inline void operator+= (const vec4 &sf) { v += sf.v; }
	inline void operator-= (const vec4 &sf) { v -= sf.v; }
	inline void operator*= (const vec4 &sf) { v *= sf.v; }
	inline void operator/= (const vec4 &sf) { v /= sf.v; }
	inline vec4 operator+ (const vec1 &sf) const { return vec4(v+vec4(sf).v); }
	inline vec4 operator- (const vec1 &sf) const { return vec4(v-vec4(sf).v); }
	inline vec4 operator* (const vec1 &sf) const { return vec4(v*vec4(sf).v); }
	inline vec4 operator/ (const vec1 &sf) const { return vec4(v/vec4(sf).v); }
	inline void operator+= (const vec1 &sf) { v += vec4(sf).v; }
	inline void operator-= (const vec1 &sf) { v -= vec4(sf).v; }
	inline void operator*= (const vec1 &sf) { v *= vec4(sf).v; }
	inline void operator/= (const vec1 &sf) { v /= vec4(sf).v; }
	inline vec4 operator^ (const vec4 &y) const { // 3D cross product
		vec4 a;
		a.s[0] = s[1]*y.s[2] - s[2]*y.s[1];
		a.s[1] = s[2]*y.s[0] - s[0]*y.s[2];
		a.s[2] = s[0]*y.s[1] - s[1]*y.s[0];
		a.s[3] = 0.0;
		return a;
	}
};

inline vec4 operator+ (const vec1 &f, const vec4 &sf) { return vec4(f)+sf; }
inline vec4 operator- (const vec1 &f, const vec4 &sf) { return vec4(f)-sf; }
inline vec4 operator* (const vec1 &f, const vec4 &sf) { return vec4(f)*sf; }
inline vec4 operator/ (const vec1 &f, const vec4 &sf) { return vec4(f)/sf; }
inline vec1 len23(vec4 x) {return (x%x);}
#endif //__VECTOR_H__
