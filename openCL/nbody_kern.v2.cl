#pragma OPENCL EXTENSION cl_khr_fp64 : enable 
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics: enable
#define TWOPI 6.28318530717958647692528676656
#ifdef PREC_DOUBLE
typedef double vec1;
typedef double2 vec2;
typedef double4 vec4;
typedef double8 vec8;
#else
typedef float vec1;
typedef float2 vec2;
typedef float4 vec4;
typedef float8 vec8;
#endif
#define CE4 (CACHE_SIZE/sizeof(vec4))
#define CE8 (CACHE_SIZE/sizeof(vec8))

typedef struct {
	uint nx;
	uint nu;
	vec1 eta;
	vec1 eps2;
	double dtn;
	double t;
} nbody_meta;

inline vec4 project(vec4 vec, vec4 nd) {
	return nd*dot(vec, nd);
}

inline vec4 plane(vec4 vec, vec4 nn) {
	return vec - project(vec, nn);
}

inline vec1 len23(vec4 v) {
	return v.x*v.x + v.y*v.y + v.z*v.z;
}

inline void block_reduce_v1dmin(const uint ti, const uint nt, __local double *tmp) {
	uint j;
	for(j = nt/2; j > 0; j >>= 1) {
		if (ti < j) { tmp[ti] = fmin(tmp[ti], tmp[ti+j]); }
		barrier(CLK_LOCAL_MEM_FENCE);
	}
}

inline void block_reduce_v4sum(const uint ti, const  uint nt, __local vec4 *tmp) {
	uint j;
	// do reduction in shared mem
	for(j = nt/2; j > 0; j >>= 1) {
		if (ti < j) { tmp[ti] += tmp[ti+j]; }
		barrier(CLK_LOCAL_MEM_FENCE);
	}
}

inline uint block_scan(const uint in, const uint ti, const uint nt, __local uint *tmp) {
	int j, ts;
	tmp[ti] = in;
	barrier(CLK_LOCAL_MEM_FENCE);
	for (j = 1; j < nt; j <<= 1) {
		ts = (ti < j) ? 0: tmp[(int)ti-j];
		barrier(CLK_LOCAL_MEM_FENCE);
		tmp[ti] += ts;
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	return tmp[ti];
}

__kernel void nbody_estep (
	__global nbody_meta *nbm,
	__write_only __global double *t1,
	__write_only __global double *t2)
{
	const uint gti = get_global_id(0), nx = nbm->nx;
	double dtn = nbm->dtn;
	if (gti < nx) {
		t1[gti] = 0.0;
		t2[gti] = dtn;
	}
}

inline vec8 nbody_accel(
	const uint gti,
	const uint nx,
	const vec1 eps,
	const vec4 xi,
	const vec4 vi,
	__read_only __global vec4 *x,
	__read_only __global vec4 *v)
{
	uint j;
	vec4 dx, dv, xj, vj;
	vec1 rp2, at;
	vec8 out;

	out = 0.0;
	for (j = 0; j < nx; j++) {
		if (j % CE8 == 0) {
			prefetch(&x[j], CE8*2);
			prefetch(&v[j], CE8*2);
		}
		if (gti != j) {
			xj = x[j];
			vj = v[j];
			dx = (vec4)(xi.xyz - xj.xyz, 0.0f);
			dv = (vec4)(vi.xyz - vj.xyz, 0.0f);
			rp2 = 1.0/(len23(dx) + eps);
			at = xj.w*rp2*sqrt(rp2);
			out.lo += at*dx;
			out.hi += at*(dv - dx*(vec1)(3.0*dot(dx,dv)*rp2));
		}
	}
	return -out;
}

__kernel void nbody_accel1(
	__constant nbody_meta *nbm,
	__read_only __global vec4 *x,
	__read_only __global vec4 *v,
	__write_only __global vec4 *a,
	__write_only __global vec4 *J)
{
	const uint gti = get_global_id(0);
	vec4 xi, vi;
	vec8 tmp;

	if (gti < nbm->nx) {
		xi = (vec4)(x[gti].xyz, 0.0f);
		vi = (vec4)(v[gti].xyz, 0.0f);
		tmp = nbody_accel(gti, nbm->nx, nbm->eps2, xi, vi, x, v);
		a[gti].xyz = tmp.s012;
		J[gti].xyz = tmp.s456;
	}
}

__kernel void nbody_f1pred(
	__constant nbody_meta *nbm,
	__read_only __global double *t1,
	__read_only __global vec4 *x,
	__read_only __global vec4 *v,
	__read_only __global vec4 *a,
	__read_only __global vec4 *J,
	__write_only __global vec4 *ox,
	__write_only __global vec4 *ov)
{
	const uint gti = get_global_id(0);
	vec4 Ji, ai, vi, xi;
	vec1 m, dt;

	if (gti < nbm->nx) {
		Ji = (vec4)(J[gti].xyz, 0.0f);
		ai = (vec4)(a[gti].xyz, 0.0f);
		vi = (vec4)(v[gti].xyz, 0.0f);
		xi = x[gti];
		dt = (vec1)(nbm->dtn - t1[gti]);
		m = xi.w;

		xi += ((Ji*dt/6.0f + ai*0.5f)*dt + vi)*dt;
		xi.w = m;
		vi += (Ji*dt*0.5f + ai)*dt;
		ox[gti] = xi;
		ov[gti].xyz = vi.xyz;
	}
}

__kernel void nbody_accel2(
	__constant nbody_meta *nbm,
	__read_only __global uint *idx,
	__global double *t1,
	__global double *t2,
	__global vec4 *x,
	__global vec4 *v,
	__global vec4 *a,
	__global vec4 *J,
	__read_only __global vec4 *x1,
	__read_only __global vec4 *v1)
{
	uint j;
	vec4 dx, dv, da, ff2, dt;
	vec4 vi, ai, Ji;
	vec4 xi1, vi1, ai1, Ji1;
	vec1 rp2, at;
	uint gti;
	const uint gid = get_global_id(0);
	double t1i, t2i, dti, dtt;
	vec8 tmp;

	if (gid < nbm->nu) {
		gti = idx[gid];
		xi1 = (vec4)(x1[gti].xyz, 0.0f);
		vi1 = (vec4)(v1[gti].xyz, 0.0f);
		tmp = nbody_accel(gti, nbm->nx, nbm->eps2, xi1, vi1, x1, v1);
		ai1 = (vec4)(tmp.lo.xyz, 0.0f);
		Ji1 = (vec4)(tmp.hi.xyz, 0.0f);

		t1i = t1[gti];
		dt = (vec1)(nbm->dtn - t1i);
		t2i = t2[gti];
		dti = t2i - t1i;
		// functionally equivalent to 
		// dx = (((Ji+Ji1)*dt/12.0 + (ai-ai1))*dt/5.0 + (vi+vi1))*dt/2.0;
		// dv = ((Ji-Ji1)*dt/6.0 + (ai+ai1))*dt/2.0;
		Ji  = (vec4)(J[gti].xyz,  0.0f);
		dx = (Ji+Ji1)*dt/12.0f;
		dv = (Ji-Ji1)*dt/6.0f;
		ai  = (vec4)(a[gti].xyz,  0.0f);
		dx = (dx + (ai-ai1))*dt*0.2f;
		dv = (dv + (ai+ai1))*dt*0.5f;
		vi  = (vec4)(v[gti].xyz,  0.0f);
		dx = (dx + (vi+vi1))*dt*0.5f;
		dv = dv + vi;

		x[gti].xyz += dx.xyz;
		v[gti].xyz = dv.xyz;
		a[gti].xyz = ai1.xyz;
		J[gti].xyz = Ji1.xyz;

		da = ai - ai1;
		ff2 = 2.0f*(3.0f*da/dt + (Ji + 2.0f*Ji1))/dt;
		dtt = sqrt(nbm->eta*sqrt(len23(ai1)/len23(ff2)));
		if (dtt < dti && dtt > 2.0*FLT_EPSILON)
			dtt = 0.5*dti;
		else if (dtt > 2.0*dti)
			dtt = 2.0*dti;
		else
			dtt = dti;
		t1[gti] = nbm->dtn;
// 		t2[gti] = tt + dtt;
		t2[gti] = nbm->dtn + pown(2.0, floor(log2(dtt)));
	}
}

__kernel void nbody_nextlist(
	__global nbody_meta *nbm,
	__read_only __global double *t2,
	__global uint *idxo,
	__local void *tmp)
{
	int cp;
	uint j, gti, cf, nx = nbm->nx;
	const uint ti = get_local_id(0);
	const uint nt = get_local_size(0);
	const uint nb = ceil((float)nx/(float)nt);
	__local int cl;
	__local double *dtmp = tmp;
	__local uint *utmp = tmp;
	double ttmp;

	ttmp = DBL_MAX;
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nx) {
			ttmp = fmin(ttmp, t2[gti]);
		}
	}
	dtmp[ti] = ttmp;
	barrier(CLK_LOCAL_MEM_FENCE);
	block_reduce_v1dmin(ti, nt, dtmp);
	ttmp = dtmp[0];

	if (ti == 0) {cl = -1;}
	for (j = 0; j < nb; j++) {
		barrier(CLK_LOCAL_MEM_FENCE);
		gti = j*nt + ti;
		cf = 0;
		if (gti < nx) {
			if (t2[gti] <= ttmp) {cf = 1;}
		}
		cp = cl + block_scan(cf, ti, nt, utmp);
		barrier(CLK_LOCAL_MEM_FENCE);
		if (cf == 1) {idxo[cp] = gti;}
		if (ti == nt-1) { cl = cp; }
	}
	if (ti == 0) {
		nbm->nu = cl;
		nbm->t = nbm->dtn;
		nbm->dtn = ttmp;
	}
}

__kernel void nbody_energy(
	__constant nbody_meta *nbm,
	__read_only __global vec4 *x,
	__read_only __global vec4 *v,
	__write_only __global vec4 *Eo)
{
	uint j;
	vec4 xi, xj;
	vec1 dx, ek, ew, m;
	vec4 E;
	const uint gti = get_global_id(0);

	if (gti < nbm->nx) {
		ek = len23(v[gti]);
		xi = x[gti];
		ew = 0.0;
		for (j = 0; j < nbm->nx; j++) {
			if (j % CE4 == 0)
				prefetch(&x[j], CE4*2);
			if (gti != j) {
				xj = x[j];
				dx = len23(xi - xj);
				ew = fma(rsqrt(dx + nbm->eps2), xj.w, ew);
			}
		}
		E = (vec4)(ek, -ew, 0.0f, 0.0f);
		Eo[gti] = xi.w*E/2.0f;
	}
}

// Find centre of mass
__kernel void array_shiftcm(
	__constant nbody_meta *nbm,
	__global vec4 *x,
	__local vec4 *tmp)
{
	uint j, gti;
	const uint ti = get_local_id(0);
	const uint nt = get_local_size(0);
	const uint nb = ceil((float)nbm->nx/(float)nt);
	vec4 ttmp, xt;
	vec1 mtmp;

	ttmp = 0.0;
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nbm->nx) {
			xt = x[gti];
			ttmp.xyz += xt.xyz*xt.w;
			ttmp.w += xt.w;
		}
	}
	tmp[ti] = ttmp;
	barrier(CLK_LOCAL_MEM_FENCE);
	block_reduce_v4sum(ti, nt, tmp);
	ttmp = tmp[0]/(tmp[0].w);
	
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nbm->nx)
			x[gti].xyz -= ttmp.xyz;
	}
}

__kernel void array_shiftbm(
	__constant nbody_meta *nbm,
	__global vec4 *v,
	__local vec4 *tmp)
{
	uint j, gti;
	const uint ti = get_local_id(0);
	const uint nt = get_local_size(0);
	const uint nb = ceil((float)nbm->nx/(float)nt);
	vec4 ttmp, xt;

	ttmp = 0.0;
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nbm->nx) {
			xt = v[gti];
			ttmp.xyz += xt.xyz;
		}
	}
	tmp[ti] = ttmp;
	barrier(CLK_LOCAL_MEM_FENCE);
	block_reduce_v4sum(ti, nt, tmp);
	ttmp = tmp[0]/((vec1)nbm->nx);
	
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nbm->nx)
			v[gti].xyz -= ttmp.xyz;
	}
}

// list reduction
// Do as much parallel reduction as possible, but reduction is serial across blocks.
// nt should be large, say 1024 for SM2 or 512 for SM1
__kernel void array_reduce_v4sum(
	const uint nx,
	__read_only __global vec4 *data,
	__write_only __global vec4 *out,
	__local vec4 *tmp)
{
	uint j, gti;
	const uint ti = get_local_id(0);
	const uint nt = get_local_size(0);
	const uint nb = ceil((float)nx/(float)nt);
	vec4 ttmp;

	ttmp = 0.0;
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nx) {
			ttmp += data[gti];
		}
	}
	tmp[ti] = ttmp;
	barrier(CLK_LOCAL_MEM_FENCE);
	block_reduce_v4sum(ti, nt, tmp);
	if (ti == 0) out[0] = tmp[0];
}

// Do as much parallel reduction as possible, but reduction is serial across blocks.
// nt should be large, say 1024.
__kernel void array_reduce_v1dmin(
	const uint nx,
	__read_only __global double *data,
	__write_only __global double *out,
	__local double *tmp)
{
	uint j, gti;
	const uint ti = get_local_id(0);
	const uint nt = get_local_size(0);
	const uint nb = ceil((float)nx/(float)nt);
	double ttmp;

	ttmp = DBL_MAX;
	for (j = 0; j < nb; j++) {
		gti = j*nt + ti;
		if (gti < nx) {
			ttmp = fmin(ttmp, data[gti]);
		}
	}
	tmp[ti] = ttmp;
	barrier(CLK_LOCAL_MEM_FENCE);
	block_reduce_v1dmin(ti, nt, tmp);
	if (ti == 0) out[0] = tmp[0];
}
