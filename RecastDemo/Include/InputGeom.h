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

#ifndef INPUTGEOM_H
#define INPUTGEOM_H

#include "ChunkyTriMesh.h"
#include "MeshLoaderObj.h"

static const int MAX_CONVEXVOL_PTS = 12;
struct ConvexVolume
{
	double verts[MAX_CONVEXVOL_PTS*3];
	double hmin, hmax;
	int nverts;
	int area;
};

struct BuildSettings
{
	// Cell size in world units
	double cellSize;
	// Cell height in world units
	double cellHeight;
	// Agent height in world units
	double agentHeight;
	// Agent radius in world units
	double agentRadius;
	// Agent max climb in world units
	double agentMaxClimb;
	// Agent max slope in degrees
	double agentMaxSlope;
	// Region minimum size in voxels.
	// regionMinSize = sqrt(regionMinArea)
	double regionMinSize;
	// Region merge size in voxels.
	// regionMergeSize = sqrt(regionMergeArea)
	double regionMergeSize;
	// Edge max length in world units
	double edgeMaxLen;
	// Edge max error in voxels
	double edgeMaxError;
	double vertsPerPoly;
	// Detail sample distance in voxels
	double detailSampleDist;
	// Detail sample max error in voxel heights.
	double detailSampleMaxError;
	// Partition type, see SamplePartitionType
	int partitionType;
	// Bounds of the area to mesh
	double navMeshBMin[3];
	double navMeshBMax[3];
	// Size of the tiles in voxels
	double tileSize;
};

class InputGeom
{
	rcChunkyTriMesh* m_chunkyMesh;
	rcMeshLoaderObj* m_mesh;
	double m_meshBMin[3], m_meshBMax[3];
	BuildSettings m_buildSettings;
	bool m_hasBuildSettings;
	
	/// @name Off-Mesh connections.
	///@{
	static const int MAX_OFFMESH_CONNECTIONS = 256;
	double m_offMeshConVerts[MAX_OFFMESH_CONNECTIONS*3*2];
	double m_offMeshConRads[MAX_OFFMESH_CONNECTIONS];
	unsigned char m_offMeshConDirs[MAX_OFFMESH_CONNECTIONS];
	unsigned char m_offMeshConAreas[MAX_OFFMESH_CONNECTIONS];
	unsigned short m_offMeshConFlags[MAX_OFFMESH_CONNECTIONS];
	unsigned int m_offMeshConId[MAX_OFFMESH_CONNECTIONS];
	int m_offMeshConCount;
	///@}

	/// @name Convex Volumes.
	///@{
	static const int MAX_VOLUMES = 256;
	ConvexVolume m_volumes[MAX_VOLUMES];
	int m_volumeCount;
	///@}
	
	bool loadMesh(class rcContext* ctx, const std::string& filepath);
	bool loadGeomSet(class rcContext* ctx, const std::string& filepath);
public:
	InputGeom();
	~InputGeom();
	
	
	bool load(class rcContext* ctx, const std::string& filepath);
	bool saveGeomSet(const BuildSettings* settings);
	
	/// Method to return static mesh data.
	const rcMeshLoaderObj* getMesh() const { return m_mesh; }
	const double* getMeshBoundsMin() const { return m_meshBMin; }
	const double* getMeshBoundsMax() const { return m_meshBMax; }
	const double* getNavMeshBoundsMin() const { return m_hasBuildSettings ? m_buildSettings.navMeshBMin : m_meshBMin; }
	const double* getNavMeshBoundsMax() const { return m_hasBuildSettings ? m_buildSettings.navMeshBMax : m_meshBMax; }
	const rcChunkyTriMesh* getChunkyMesh() const { return m_chunkyMesh; }
	const BuildSettings* getBuildSettings() const { return m_hasBuildSettings ? &m_buildSettings : 0; }
	bool raycastMesh(double* src, double* dst, double& tmin);

	/// @name Off-Mesh connections.
	///@{
	int getOffMeshConnectionCount() const { return m_offMeshConCount; }
	const double* getOffMeshConnectionVerts() const { return m_offMeshConVerts; }
	const double* getOffMeshConnectionRads() const { return m_offMeshConRads; }
	const unsigned char* getOffMeshConnectionDirs() const { return m_offMeshConDirs; }
	const unsigned char* getOffMeshConnectionAreas() const { return m_offMeshConAreas; }
	const unsigned short* getOffMeshConnectionFlags() const { return m_offMeshConFlags; }
	const unsigned int* getOffMeshConnectionId() const { return m_offMeshConId; }
	void addOffMeshConnection(const double* spos, const double* epos, const double rad,
							  unsigned char bidir, unsigned char area, unsigned short flags);
	void deleteOffMeshConnection(int i);
	void drawOffMeshConnections(struct duDebugDraw* dd, bool hilight = false);
	///@}

	/// @name Box Volumes.
	///@{
	int getConvexVolumeCount() const { return m_volumeCount; }
	const ConvexVolume* getConvexVolumes() const { return m_volumes; }
	void addConvexVolume(const double* verts, const int nverts,
						 const double minh, const double maxh, unsigned char area);
	void deleteConvexVolume(int i);
	void drawConvexVolumes(struct duDebugDraw* dd, bool hilight = false);
	///@}
	
private:
	// Explicitly disabled copy constructor and copy assignment operator.
	InputGeom(const InputGeom&);
	InputGeom& operator=(const InputGeom&);
};

#endif // INPUTGEOM_H
