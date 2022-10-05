//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#include "DetourObstacleAvoidance.h"
#include "DetourCommon.h"
#include "DetourMath.h"
#include "DetourAlloc.h"
#include "DetourAssert.h"
#include <string.h>
#include <float.h>
#include <new>

static const double DT_PI = 3.141592653589793238462643383279502884;

static int sweepCircleCircle(const double* c0, const double r0, const double* v,
							 const double* c1, const double r1,
							 double& tmin, double& tmax)
{
	static const double EPS = 0.0001;
	double s[3];
	dtVsub(s,c1,c0);
	double r = r0+r1;
	double c = dtVdot2D(s,s) - r*r;
	double a = dtVdot2D(v,v);
	if (a < EPS) return 0;	// not moving
	
	// Overlap, calc time to exit.
	double b = dtVdot2D(v,s);
	double d = b*b - a*c;
	if (d < 0.0) return 0; // no intersection.
	a = 1.0 / a;
	const double rd = dtMathSqrtf(d);
	tmin = (b - rd) * a;
	tmax = (b + rd) * a;
	return 1;
}

static int isectRaySeg(const double* ap, const double* u,
					   const double* bp, const double* bq,
					   double& t)
{
	double v[3], w[3];
	dtVsub(v,bq,bp);
	dtVsub(w,ap,bp);
	double d = dtVperp2D(u,v);
	if (dtMathFabsf(d) < 1e-6) return 0;
	d = 1.0/d;
	t = dtVperp2D(v,w) * d;
	if (t < 0 || t > 1) return 0;
	double s = dtVperp2D(u,w) * d;
	if (s < 0 || s > 1) return 0;
	return 1;
}



dtObstacleAvoidanceDebugData* dtAllocObstacleAvoidanceDebugData()
{
	void* mem = dtAlloc(sizeof(dtObstacleAvoidanceDebugData), DT_ALLOC_PERM);
	if (!mem) return 0;
	return new(mem) dtObstacleAvoidanceDebugData;
}

void dtFreeObstacleAvoidanceDebugData(dtObstacleAvoidanceDebugData* ptr)
{
	if (!ptr) return;
	ptr->~dtObstacleAvoidanceDebugData();
	dtFree(ptr);
}


dtObstacleAvoidanceDebugData::dtObstacleAvoidanceDebugData() :
	m_nsamples(0),
	m_maxSamples(0),
	m_vel(0),
	m_ssize(0),
	m_pen(0),
	m_vpen(0),
	m_vcpen(0),
	m_spen(0),
	m_tpen(0)
{
}

dtObstacleAvoidanceDebugData::~dtObstacleAvoidanceDebugData()
{
	dtFree(m_vel);
	dtFree(m_ssize);
	dtFree(m_pen);
	dtFree(m_vpen);
	dtFree(m_vcpen);
	dtFree(m_spen);
	dtFree(m_tpen);
}
		
bool dtObstacleAvoidanceDebugData::init(const int maxSamples)
{
	dtAssert(maxSamples);
	m_maxSamples = maxSamples;

	m_vel = (double*)dtAlloc(sizeof(double)*3*m_maxSamples, DT_ALLOC_PERM);
	if (!m_vel)
		return false;
	m_pen = (double*)dtAlloc(sizeof(double)*m_maxSamples, DT_ALLOC_PERM);
	if (!m_pen)
		return false;
	m_ssize = (double*)dtAlloc(sizeof(double)*m_maxSamples, DT_ALLOC_PERM);
	if (!m_ssize)
		return false;
	m_vpen = (double*)dtAlloc(sizeof(double)*m_maxSamples, DT_ALLOC_PERM);
	if (!m_vpen)
		return false;
	m_vcpen = (double*)dtAlloc(sizeof(double)*m_maxSamples, DT_ALLOC_PERM);
	if (!m_vcpen)
		return false;
	m_spen = (double*)dtAlloc(sizeof(double)*m_maxSamples, DT_ALLOC_PERM);
	if (!m_spen)
		return false;
	m_tpen = (double*)dtAlloc(sizeof(double)*m_maxSamples, DT_ALLOC_PERM);
	if (!m_tpen)
		return false;
	
	return true;
}

void dtObstacleAvoidanceDebugData::reset()
{
	m_nsamples = 0;
}

void dtObstacleAvoidanceDebugData::addSample(const double* vel, const double ssize, const double pen,
											 const double vpen, const double vcpen, const double spen, const double tpen)
{
	if (m_nsamples >= m_maxSamples)
		return;
	dtAssert(m_vel);
	dtAssert(m_ssize);
	dtAssert(m_pen);
	dtAssert(m_vpen);
	dtAssert(m_vcpen);
	dtAssert(m_spen);
	dtAssert(m_tpen);
	dtVcopy(&m_vel[m_nsamples*3], vel);
	m_ssize[m_nsamples] = ssize;
	m_pen[m_nsamples] = pen;
	m_vpen[m_nsamples] = vpen;
	m_vcpen[m_nsamples] = vcpen;
	m_spen[m_nsamples] = spen;
	m_tpen[m_nsamples] = tpen;
	m_nsamples++;
}

static void normalizeArray(double* arr, const int n)
{
	// Normalize penaly range.
	double minPen = DBL_MAX;
	double maxPen = -DBL_MAX;
	for (int i = 0; i < n; ++i)
	{
		minPen = dtMin(minPen, arr[i]);
		maxPen = dtMax(maxPen, arr[i]);
	}
	const double penRange = maxPen-minPen;
	const double s = penRange > 0.001 ? (1.0 / penRange) : 1;
	for (int i = 0; i < n; ++i)
		arr[i] = dtClamp((arr[i]-minPen)*s, 0.0, 1.0);
}

void dtObstacleAvoidanceDebugData::normalizeSamples()
{
	normalizeArray(m_pen, m_nsamples);
	normalizeArray(m_vpen, m_nsamples);
	normalizeArray(m_vcpen, m_nsamples);
	normalizeArray(m_spen, m_nsamples);
	normalizeArray(m_tpen, m_nsamples);
}


dtObstacleAvoidanceQuery* dtAllocObstacleAvoidanceQuery()
{
	void* mem = dtAlloc(sizeof(dtObstacleAvoidanceQuery), DT_ALLOC_PERM);
	if (!mem) return 0;
	return new(mem) dtObstacleAvoidanceQuery;
}

void dtFreeObstacleAvoidanceQuery(dtObstacleAvoidanceQuery* ptr)
{
	if (!ptr) return;
	ptr->~dtObstacleAvoidanceQuery();
	dtFree(ptr);
}


dtObstacleAvoidanceQuery::dtObstacleAvoidanceQuery() :
	m_invHorizTime(0),
	m_vmax(0),
	m_invVmax(0),
	m_maxCircles(0),
	m_circles(0),
	m_ncircles(0),
	m_maxSegments(0),
	m_segments(0),
	m_nsegments(0)
{
}

dtObstacleAvoidanceQuery::~dtObstacleAvoidanceQuery()
{
	dtFree(m_circles);
	dtFree(m_segments);
}

bool dtObstacleAvoidanceQuery::init(const int maxCircles, const int maxSegments)
{
	m_maxCircles = maxCircles;
	m_ncircles = 0;
	m_circles = (dtObstacleCircle*)dtAlloc(sizeof(dtObstacleCircle)*m_maxCircles, DT_ALLOC_PERM);
	if (!m_circles)
		return false;
	memset(m_circles, 0, sizeof(dtObstacleCircle)*m_maxCircles);

	m_maxSegments = maxSegments;
	m_nsegments = 0;
	m_segments = (dtObstacleSegment*)dtAlloc(sizeof(dtObstacleSegment)*m_maxSegments, DT_ALLOC_PERM);
	if (!m_segments)
		return false;
	memset(m_segments, 0, sizeof(dtObstacleSegment)*m_maxSegments);
	
	return true;
}

void dtObstacleAvoidanceQuery::reset()
{
	m_ncircles = 0;
	m_nsegments = 0;
}

void dtObstacleAvoidanceQuery::addCircle(const double* pos, const double rad,
										 const double* vel, const double* dvel)
{
	if (m_ncircles >= m_maxCircles)
		return;
		
	dtObstacleCircle* cir = &m_circles[m_ncircles++];
	dtVcopy(cir->p, pos);
	cir->rad = rad;
	dtVcopy(cir->vel, vel);
	dtVcopy(cir->dvel, dvel);
}

void dtObstacleAvoidanceQuery::addSegment(const double* p, const double* q)
{
	if (m_nsegments >= m_maxSegments)
		return;
	
	dtObstacleSegment* seg = &m_segments[m_nsegments++];
	dtVcopy(seg->p, p);
	dtVcopy(seg->q, q);
}

void dtObstacleAvoidanceQuery::prepare(const double* pos, const double* dvel)
{
	// Prepare obstacles
	for (int i = 0; i < m_ncircles; ++i)
	{
		dtObstacleCircle* cir = &m_circles[i];
		
		// Side
		const double* pa = pos;
		const double* pb = cir->p;
		
		const double orig[3] = {0,0,0};
		double dv[3];
		dtVsub(cir->dp,pb,pa);
		dtVnormalize(cir->dp);
		dtVsub(dv, cir->dvel, dvel);
		
		const double a = dtTriArea2D(orig, cir->dp,dv);
		if (a < 0.01)
		{
			cir->np[0] = -cir->dp[2];
			cir->np[2] = cir->dp[0];
		}
		else
		{
			cir->np[0] = cir->dp[2];
			cir->np[2] = -cir->dp[0];
		}
	}	

	for (int i = 0; i < m_nsegments; ++i)
	{
		dtObstacleSegment* seg = &m_segments[i];
		
		// Precalc if the agent is really close to the segment.
		const double r = 0.01;
		double t;
		seg->touch = dtDistancePtSegSqr2D(pos, seg->p, seg->q, t) < dtSqr(r);
	}	
}


/* Calculate the collision penalty for a given velocity vector
 * 
 * @param vcand sampled velocity
 * @param dvel desired velocity
 * @param minPenalty threshold penalty for early out
 */
double dtObstacleAvoidanceQuery::processSample(const double* vcand, const double cs,
											  const double* pos, const double rad,
											  const double* vel, const double* dvel,
											  const double minPenalty,
											  dtObstacleAvoidanceDebugData* debug)
{
	// penalty for straying away from the desired and current velocities
	const double vpen = m_params.weightDesVel * (dtVdist2D(vcand, dvel) * m_invVmax);
	const double vcpen = m_params.weightCurVel * (dtVdist2D(vcand, vel) * m_invVmax);

	// find the threshold hit time to bail out based on the early out penalty
	// (see how the penalty is calculated below to understnad)
	double minPen = minPenalty - vpen - vcpen;
	double tThresold = (m_params.weightToi / minPen - 0.1) * m_params.horizTime;
	if (tThresold - m_params.horizTime > -FLT_EPSILON)
		return minPenalty; // already too much

	// Find min time of impact and exit amongst all obstacles.
	double tmin = m_params.horizTime;
	double side = 0;
	int nside = 0;
	const double zero = 0.0;
    const double one = 1.0;
	for (int i = 0; i < m_ncircles; ++i)
	{
		const dtObstacleCircle* cir = &m_circles[i];
			
		// RVO
		double vab[3];
		dtVscale(vab, vcand, 2);
		dtVsub(vab, vab, vel);
		dtVsub(vab, vab, cir->vel);
		
		// Side
		side += dtClamp(dtMin(dtVdot2D(cir->dp,vab)*0.5+0.5, dtVdot2D(cir->np,vab)*2), zero, one);
		nside++;
		
		double htmin = 0, htmax = 0;
		if (!sweepCircleCircle(pos,rad, vab, cir->p,cir->rad, htmin, htmax))
			continue;
		
		// Handle overlapping obstacles.
		if (htmin < 0.0 && htmax > 0.0)
		{
			// Avoid more when overlapped.
			htmin = -htmin * 0.5;
		}
		
		if (htmin >= 0.0)
		{
			// The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
			if (htmin < tmin)
			{
				tmin = htmin;
				if (tmin < tThresold)
					return minPenalty;
			}
		}
	}

	for (int i = 0; i < m_nsegments; ++i)
	{
		const dtObstacleSegment* seg = &m_segments[i];
		double htmin = 0;
		
		if (seg->touch)
		{
			// Special case when the agent is very close to the segment.
			double sdir[3], snorm[3];
			dtVsub(sdir, seg->q, seg->p);
			snorm[0] = -sdir[2];
			snorm[2] = sdir[0];
			// If the velocity is pointing towards the segment, no collision.
			if (dtVdot2D(snorm, vcand) < 0.0)
				continue;
			// Else immediate collision.
			htmin = 0.0;
		}
		else
		{
			if (!isectRaySeg(pos, vcand, seg->p, seg->q, htmin))
				continue;
		}
		
		// Avoid less when facing walls.
		htmin *= 2.0;
		
		// The closest obstacle is somewhere ahead of us, keep track of nearest obstacle.
		if (htmin < tmin)
		{
			tmin = htmin;
			if (tmin < tThresold)
				return minPenalty;
		}
	}
	
	// Normalize side bias, to prevent it dominating too much.
	if (nside)
		side /= nside;
	
	const double spen = m_params.weightSide * side;
	const double tpen = m_params.weightToi * (1.0/(0.1+tmin*m_invHorizTime));
	
	const double penalty = vpen + vcpen + spen + tpen;
	
	// Store different penalties for debug viewing
	if (debug)
		debug->addSample(vcand, cs, penalty, vpen, vcpen, spen, tpen);
	
	return penalty;
}

int dtObstacleAvoidanceQuery::sampleVelocityGrid(const double* pos, const double rad, const double vmax,
												 const double* vel, const double* dvel, double* nvel,
												 const dtObstacleAvoidanceParams* params,
												 dtObstacleAvoidanceDebugData* debug)
{
	prepare(pos, dvel);
	
	memcpy(&m_params, params, sizeof(dtObstacleAvoidanceParams));
	m_invHorizTime = 1.0 / m_params.horizTime;
	m_vmax = vmax;
	m_invVmax = vmax > 0 ? 1.0 / vmax : DBL_MAX;
	
	dtVset(nvel, 0,0,0);
	
	if (debug)
		debug->reset();

	const double cvx = dvel[0] * m_params.velBias;
	const double cvz = dvel[2] * m_params.velBias;
	const double cs = vmax * 2 * (1 - m_params.velBias) / (double)(m_params.gridSize-1);
	const double half = (m_params.gridSize-1)*cs*0.5;
		
	double minPenalty = DBL_MAX;
	int ns = 0;
		
	for (int y = 0; y < m_params.gridSize; ++y)
	{
		for (int x = 0; x < m_params.gridSize; ++x)
		{
			double vcand[3];
			vcand[0] = cvx + x*cs - half;
			vcand[1] = 0;
			vcand[2] = cvz + y*cs - half;
			
			if (dtSqr(vcand[0])+dtSqr(vcand[2]) > dtSqr(vmax+cs/2)) continue;
			
			const double penalty = processSample(vcand, cs, pos,rad,vel,dvel, minPenalty, debug);
			ns++;
			if (penalty < minPenalty)
			{
				minPenalty = penalty;
				dtVcopy(nvel, vcand);
			}
		}
	}
	
	return ns;
}


// vector normalization that ignores the y-component.
inline void dtNormalize2D(double* v)
{
	double d = dtMathSqrtf(v[0] * v[0] + v[2] * v[2]);
	if (d==0)
		return;
	d = 1.0 / d;
	v[0] *= d;
	v[2] *= d;
}

// vector normalization that ignores the y-component.
inline void dtRorate2D(double* dest, const double* v, double ang)
{
	double c = cos(ang);
	double s = sin(ang);
	dest[0] = v[0]*c - v[2]*s;
	dest[2] = v[0]*s + v[2]*c;
	dest[1] = v[1];
}


int dtObstacleAvoidanceQuery::sampleVelocityAdaptive(const double* pos, const double rad, const double vmax,
													 const double* vel, const double* dvel, double* nvel,
													 const dtObstacleAvoidanceParams* params,
													 dtObstacleAvoidanceDebugData* debug)
{
	prepare(pos, dvel);
	
	memcpy(&m_params, params, sizeof(dtObstacleAvoidanceParams));
	m_invHorizTime = 1.0 / m_params.horizTime;
	m_vmax = vmax;
	m_invVmax = vmax > 0 ? 1.0 / vmax : DBL_MAX;
	
	dtVset(nvel, 0,0,0);
	
	if (debug)
		debug->reset();

	// Build sampling pattern aligned to desired velocity.
	double pat[(DT_MAX_PATTERN_DIVS*DT_MAX_PATTERN_RINGS+1)*2];
	int npat = 0;

	const int ndivs = (int)m_params.adaptiveDivs;
	const int nrings= (int)m_params.adaptiveRings;
	const int depth = (int)m_params.adaptiveDepth;
	
	const int nd = dtClamp(ndivs, 1, DT_MAX_PATTERN_DIVS);
	const int nr = dtClamp(nrings, 1, DT_MAX_PATTERN_RINGS);
	const double da = (1.0/nd) * DT_PI*2;
	const double ca = cos(da);
	const double sa = sin(da);

	// desired direction
	double ddir[6];
	dtVcopy(ddir, dvel);
	dtNormalize2D(ddir);
	dtRorate2D (ddir+3, ddir, da*0.5); // rotated by da/2

	// Always add sample at zero
	pat[npat*2+0] = 0;
	pat[npat*2+1] = 0;
	npat++;
	
	for (int j = 0; j < nr; ++j)
	{
		const double r = (double)(nr-j)/(double)nr;
		pat[npat*2+0] = ddir[(j%2)*3] * r;
		pat[npat*2+1] = ddir[(j%2)*3+2] * r;
		double* last1 = pat + npat*2;
		double* last2 = last1;
		npat++;

		for (int i = 1; i < nd-1; i+=2)
		{
			// get next point on the "right" (rotate CW)
			pat[npat*2+0] = last1[0]*ca + last1[1]*sa;
			pat[npat*2+1] = -last1[0]*sa + last1[1]*ca;
			// get next point on the "left" (rotate CCW)
			pat[npat*2+2] = last2[0]*ca - last2[1]*sa;
			pat[npat*2+3] = last2[0]*sa + last2[1]*ca;

			last1 = pat + npat*2;
			last2 = last1 + 2;
			npat += 2;
		}

		if ((nd&1) == 0)
		{
			pat[npat*2+2] = last2[0]*ca - last2[1]*sa;
			pat[npat*2+3] = last2[0]*sa + last2[1]*ca;
			npat++;
		}
	}


	// Start sampling.
	double cr = vmax * (1.0 - m_params.velBias);
	double res[3];
	dtVset(res, dvel[0] * m_params.velBias, 0, dvel[2] * m_params.velBias);
	int ns = 0;

	for (int k = 0; k < depth; ++k)
	{
		double minPenalty = DBL_MAX;
		double bvel[3];
		dtVset(bvel, 0,0,0);
		
		for (int i = 0; i < npat; ++i)
		{
			double vcand[3];
			vcand[0] = res[0] + pat[i*2+0]*cr;
			vcand[1] = 0;
			vcand[2] = res[2] + pat[i*2+1]*cr;
			
			if (dtSqr(vcand[0])+dtSqr(vcand[2]) > dtSqr(vmax+0.001)) continue;
			
			const double penalty = processSample(vcand,cr/10, pos,rad,vel,dvel, minPenalty, debug);
			ns++;
			if (penalty < minPenalty)
			{
				minPenalty = penalty;
				dtVcopy(bvel, vcand);
			}
		}

		dtVcopy(res, bvel);

		cr *= 0.5;
	}	
	
	dtVcopy(nvel, res);
	
	return ns;
}
