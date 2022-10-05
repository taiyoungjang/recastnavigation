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

#ifndef DETOUROBSTACLEAVOIDANCE_H
#define DETOUROBSTACLEAVOIDANCE_H

struct dtObstacleCircle
{
	double p[3];				///< Position of the obstacle
	double vel[3];			///< Velocity of the obstacle
	double dvel[3];			///< Velocity of the obstacle
	double rad;				///< Radius of the obstacle
	double dp[3], np[3];		///< Use for side selection during sampling.
};

struct dtObstacleSegment
{
	double p[3], q[3];		///< End points of the obstacle segment
	bool touch;
};


class dtObstacleAvoidanceDebugData
{
public:
	dtObstacleAvoidanceDebugData();
	~dtObstacleAvoidanceDebugData();
	
	bool init(const int maxSamples);
	void reset();
	void addSample(const double* vel, const double ssize, const double pen,
				   const double vpen, const double vcpen, const double spen, const double tpen);
	
	void normalizeSamples();
	
	inline int getSampleCount() const { return m_nsamples; }
	inline const double* getSampleVelocity(const int i) const { return &m_vel[i*3]; }
	inline double getSampleSize(const int i) const { return m_ssize[i]; }
	inline double getSamplePenalty(const int i) const { return m_pen[i]; }
	inline double getSampleDesiredVelocityPenalty(const int i) const { return m_vpen[i]; }
	inline double getSampleCurrentVelocityPenalty(const int i) const { return m_vcpen[i]; }
	inline double getSamplePreferredSidePenalty(const int i) const { return m_spen[i]; }
	inline double getSampleCollisionTimePenalty(const int i) const { return m_tpen[i]; }

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	dtObstacleAvoidanceDebugData(const dtObstacleAvoidanceDebugData&);
	dtObstacleAvoidanceDebugData& operator=(const dtObstacleAvoidanceDebugData&);

	int m_nsamples;
	int m_maxSamples;
	double* m_vel;
	double* m_ssize;
	double* m_pen;
	double* m_vpen;
	double* m_vcpen;
	double* m_spen;
	double* m_tpen;
};

dtObstacleAvoidanceDebugData* dtAllocObstacleAvoidanceDebugData();
void dtFreeObstacleAvoidanceDebugData(dtObstacleAvoidanceDebugData* ptr);


static const int DT_MAX_PATTERN_DIVS = 32;	///< Max numver of adaptive divs.
static const int DT_MAX_PATTERN_RINGS = 4;	///< Max number of adaptive rings.

struct dtObstacleAvoidanceParams
{
	double velBias;
	double weightDesVel;
	double weightCurVel;
	double weightSide;
	double weightToi;
	double horizTime;
	unsigned char gridSize;	///< grid
	unsigned char adaptiveDivs;	///< adaptive
	unsigned char adaptiveRings;	///< adaptive
	unsigned char adaptiveDepth;	///< adaptive
};

class dtObstacleAvoidanceQuery
{
public:
	dtObstacleAvoidanceQuery();
	~dtObstacleAvoidanceQuery();
	
	bool init(const int maxCircles, const int maxSegments);
	
	void reset();

	void addCircle(const double* pos, const double rad,
				   const double* vel, const double* dvel);
				   
	void addSegment(const double* p, const double* q);

	int sampleVelocityGrid(const double* pos, const double rad, const double vmax,
						   const double* vel, const double* dvel, double* nvel,
						   const dtObstacleAvoidanceParams* params,
						   dtObstacleAvoidanceDebugData* debug = 0);

	int sampleVelocityAdaptive(const double* pos, const double rad, const double vmax,
							   const double* vel, const double* dvel, double* nvel,
							   const dtObstacleAvoidanceParams* params, 
							   dtObstacleAvoidanceDebugData* debug = 0);
	
	inline int getObstacleCircleCount() const { return m_ncircles; }
	const dtObstacleCircle* getObstacleCircle(const int i) { return &m_circles[i]; }

	inline int getObstacleSegmentCount() const { return m_nsegments; }
	const dtObstacleSegment* getObstacleSegment(const int i) { return &m_segments[i]; }

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	dtObstacleAvoidanceQuery(const dtObstacleAvoidanceQuery&);
	dtObstacleAvoidanceQuery& operator=(const dtObstacleAvoidanceQuery&);

	void prepare(const double* pos, const double* dvel);

	double processSample(const double* vcand, const double cs,
						const double* pos, const double rad,
						const double* vel, const double* dvel,
						const double minPenalty,
						dtObstacleAvoidanceDebugData* debug);

	dtObstacleAvoidanceParams m_params;
	double m_invHorizTime;
	double m_vmax;
	double m_invVmax;

	int m_maxCircles;
	dtObstacleCircle* m_circles;
	int m_ncircles;

	int m_maxSegments;
	dtObstacleSegment* m_segments;
	int m_nsegments;
};

dtObstacleAvoidanceQuery* dtAllocObstacleAvoidanceQuery();
void dtFreeObstacleAvoidanceQuery(dtObstacleAvoidanceQuery* ptr);


#endif // DETOUROBSTACLEAVOIDANCE_H
